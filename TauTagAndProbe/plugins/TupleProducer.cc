/*! Creates tuple for tau analysis.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#include "Compression.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TauTriggerTools/Common/interface/CutTools.h"
#include "TauTriggerTools/Common/interface/GenTruthTools.h"
#include "TauTriggerTools/Common/interface/PatHelpers.h"
#include "TauTriggerTools/Common/interface/TriggerDescriptor.h"
#include "TauTriggerTools/TauTagAndProbe/interface/SummaryTuple.h"

namespace tau_trigger {

struct TupleProducerData {
    using Mutex = EventTuple::Mutex;
    using LockGuard = std::lock_guard<Mutex>;
    using SelectionHist = root_ext::SmartHistogram<cuts::ObjectSelector>;

    std::unique_ptr<EventTuple> eventTuple;
    std::unique_ptr<SelectionHist> selection;

    TupleProducerData(TFile& file)
    {
        eventTuple = std::make_unique<EventTuple>("events", &file, false);
        selection = std::make_unique<SelectionHist>("producer_selection");
        selection->SetOutputDirectory(&file);
    }
};

class TupleProducer : public edm::stream::EDProducer<edm::GlobalCache<TupleProducerData>> {
public:
    using SelectionHist = TupleProducerData::SelectionHist;
    using Cutter = cuts::Cutter<>;
    using exception = analysis::exception;

    TupleProducer(const edm::ParameterSet& cfg, const TupleProducerData* producerData) :
        btagThreshold(cfg.getParameter<double>("btagThreshold")),
        isMC(cfg.getParameter<bool>("isMC")),
        triggerProcess(cfg.getParameter<std::string>("triggerProcess")),
        genEvent_token(consumeIT<GenEventInfoProduct>(cfg, "genEvent", false)),
        genParticles_token(consumeIT<reco::GenParticleCollection>(cfg, "genParticles", false)),
        puInfo_token(consumeIT<std::vector<PileupSummaryInfo>>(cfg, "puInfo", false)),
        vertices_token(consumeIT<std::vector<reco::Vertex>>(cfg, "vertices")),
        taus_token(consumeIT<pat::TauCollection>(cfg, "taus")),
        jets_token(consumeIT<pat::JetCollection>(cfg, "jets")),
        met_token(consumeIT<pat::METCollection>(cfg, "met")),
        triggerResults_token(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", triggerProcess))),
        triggerObjects_token(consumeIT<pat::TriggerObjectStandAloneCollection>(cfg, "triggerObjects")),
        l1Taus_token(consumeIT<l1t::TauBxCollection>(cfg, "l1Taus")),
        triggerDescriptors(cfg.getParameter<edm::VParameterSet>("hltPaths")),
        data(*producerData)
    {
        produces<bool>();
    }

    static std::unique_ptr<TupleProducerData> initializeGlobalCache(const edm::ParameterSet&)
    {
        TFile& file = edm::Service<TFileService>()->file();
        file.SetCompressionAlgorithm(ROOT::kLZ4);
        file.SetCompressionLevel(4);
        return std::make_unique<TupleProducerData>(file);
    }

    static void globalEndJob(TupleProducerData* data)
    {
        TupleProducerData::LockGuard lock(data->eventTuple->GetMutex());
        data->eventTuple->Write();
        data->selection->WriteRootObject();
    }

private:
    static constexpr float default_value = ::tau_trigger::DefaultFillValue<float>();
    static constexpr int default_int_value = ::tau_trigger::DefaultFillValue<int>();
    static constexpr double deltaR2Thr = 0.5*0.5;

    template<typename T>
    edm::EDGetTokenT<T> consumeIT(const edm::ParameterSet& cfg, const std::string& name, bool always = true)
    {
        if(always)
            return consumes<T>(cfg.getParameter<edm::InputTag>(name));
        return mayConsume<T>(cfg.getParameter<edm::InputTag>(name));
    }

    virtual void beginRun(const edm::Run& run, const edm::EventSetup& setup)
    {
        HLTConfigProvider hltConfigProvider;
        bool changedConfig;
        if(!hltConfigProvider.init(run, setup, triggerProcess, changedConfig))
            throw exception("Unable to initialize HLTConfigProvider.");
        triggerDescriptors.updateGlobalIndices(hltConfigProvider.triggerNames());
    }

    virtual void produce(edm::Event& event, const edm::EventSetup&) override
    {
        event.put(std::make_unique<bool>(true));
        TupleProducerData::LockGuard lock(data.eventTuple->GetMutex());
        try {
            Cutter cut(data.selection.get());
            fillTuple(event, cut);
        } catch(cuts::cut_failed&) {}
        data.selection->fill_selection();
    }

    void fillTuple(edm::Event& event, Cutter& cut)
    {
        EventTuple& eventTuple = *data.eventTuple;
        // std::cout << "start" << std::endl;

        cut(true, "total");
        eventTuple().run  = event.id().run();
        eventTuple().lumi = event.id().luminosityBlock();
        eventTuple().evt  = event.id().event();

        edm::Handle<std::vector<reco::Vertex>> vertices;
        event.getByToken(vertices_token, vertices);
        eventTuple().npv = static_cast<int>(vertices->size());
        // std::cout << "tokens" << std::endl;

        edm::Handle<std::vector<reco::GenParticle>> hGenParticles;

        if(isMC) {
            edm::Handle<GenEventInfoProduct> genEvent;
            event.getByToken(genEvent_token, genEvent);
            eventTuple().genEventWeight = static_cast<float>(genEvent->weight());

            edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
            event.getByToken(puInfo_token, puInfo);
            eventTuple().npu = gen_truth::GetNumberOfPileUpInteractions(puInfo);

            event.getByToken(genParticles_token, hGenParticles);
        }

        // std::cout << "isMC passes" << std::endl;


        auto genParticles = hGenParticles.isValid() ? hGenParticles.product() : nullptr;
        std::vector<gen_truth::LeptonMatchResult> genLeptons;
        if(genParticles)
            genLeptons = gen_truth::CollectGenLeptons(*genParticles);

        edm::Handle<edm::TriggerResults> triggerResults;
        event.getByToken(triggerResults_token, triggerResults);
        const edm::TriggerNames& triggerNames = event.triggerNames(*triggerResults);
        edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
        event.getByToken(triggerObjects_token, triggerObjects);
        edm::Handle<l1t::TauBxCollection> l1Taus;
        event.getByToken(l1Taus_token, l1Taus);

        // std::cout << "more tokens" << std::endl;


        edm::Handle<pat::METCollection> metCollection;
        event.getByToken(met_token, metCollection);
        const pat::MET& met = metCollection->at(0);
        const LorentzVectorM met_p4(met.pt(), 0, met.phi(), 0);
        eventTuple().met_pt = static_cast<float>(met.pt());
        eventTuple().met_phi = static_cast<float>(met.phi());

        // std::cout << "met filled" << std::endl;


        edm::Handle<pat::TauCollection> taus;
        event.getByToken(taus_token, taus);

        // std::cout << "tau token" << std::endl;


        edm::Handle<pat::JetCollection> jets;
        event.getByToken(jets_token, jets);

        // std::cout << "jet token" << std::endl;


        const auto& selected_tau_pairs = CollectTauPairs(*taus, genLeptons, deltaR2Thr, cut);
        // std::cout << "cutting on pairs" << std::endl;
        cut(!selected_tau_pairs.empty(), "has_tau");
        // std::cout << "cutting on pairs successful" << std::endl;

        bool has_good_tau_pair = false;        

        for(const auto& tau_pair_entry : selected_tau_pairs) {
            
            // bool bveto = false;
            std::map<size_t, TriggerObjectMatchResult> trig_obj;
            // for(const auto& tau_entry : tau_pair_entry){
            //     const pat::Tau* tau = tau_entry.reco_tau;
            //     const auto& gen_tau = tau_entry.gen_tau;
            //     const LorentzVectorM tau_ref_p4 = tau ? tau->polarP4() : LorentzVectorM(gen_tau.visible_p4);
            //     if(btagThreshold > 0  && !PassBtagVeto(tau_ref_p4, *jets, btagThreshold, deltaR2Thr)){
            //         bveto = true;
            //         break;
            //     }
            // }
            // if(bveto) continue;
            for(const auto& tau_entry : tau_pair_entry){
                // std::cout << "new tau" << std::endl;
                const pat::Tau* tau = tau_entry.reco_tau;
                const auto& gen_tau = tau_entry.gen_tau;
                const bool has_gen_tau = gen_tau.match != GenLeptonMatch::NoMatch;
                const LorentzVectorM tau_ref_p4 = tau ? tau->polarP4() : LorentzVectorM(gen_tau.visible_p4);
                if(!tau && !has_gen_tau)
                    throw exception("Inconsistent tau entry");

                eventTuple().tau_sel.push_back(tau_entry.selection);
                eventTuple().tau_pt.push_back(tau ? static_cast<float>(tau->polarP4().pt()) : default_value);
                eventTuple().tau_eta.push_back(tau ? static_cast<float>(tau->polarP4().eta()) : default_value);
                eventTuple().tau_phi.push_back(tau ? static_cast<float>(tau->polarP4().phi()) : default_value);
                eventTuple().tau_mass.push_back(tau ? static_cast<float>(tau->polarP4().mass()) : default_value);
                eventTuple().tau_charge.push_back(tau ? tau->charge() : default_int_value);

                eventTuple().tau_gen_match.push_back(static_cast<int>(gen_tau.match));
                eventTuple().tau_gen_charge.push_back(has_gen_tau ? gen_tau.gen_particle_firstCopy->charge() : default_int_value);
                eventTuple().tau_gen_vis_pt.push_back(has_gen_tau ? static_cast<float>(gen_tau.visible_p4.pt()) : default_value);
                eventTuple().tau_gen_vis_eta.push_back(has_gen_tau ? static_cast<float>(gen_tau.visible_p4.eta()) : default_value);
                eventTuple().tau_gen_vis_phi.push_back(has_gen_tau ? static_cast<float>(gen_tau.visible_p4.phi()) : default_value);
                eventTuple().tau_gen_vis_mass.push_back(has_gen_tau ? static_cast<float>(gen_tau.visible_p4.mass()) : default_value);
                eventTuple().tau_gen_rad_pt.push_back(has_gen_tau ? static_cast<float>(gen_tau.visible_rad_p4.pt()) : default_value);
                eventTuple().tau_gen_rad_eta.push_back(has_gen_tau ? static_cast<float>(gen_tau.visible_rad_p4.eta())
                                                        : default_value);
                eventTuple().tau_gen_rad_phi.push_back(has_gen_tau ? static_cast<float>(gen_tau.visible_rad_p4.phi())
                                                        : default_value);
                eventTuple().tau_gen_rad_energy.push_back(has_gen_tau ? static_cast<float>(gen_tau.visible_rad_p4.energy())
                                                            : default_value);
                eventTuple().tau_gen_n_charged_hadrons.push_back(has_gen_tau ? static_cast<int>(gen_tau.n_charged_hadrons)
                                                                    : default_int_value);
                eventTuple().tau_gen_n_neutral_hadrons.push_back(has_gen_tau ? static_cast<int>(gen_tau.n_neutral_hadrons)
                                                                    : default_int_value);
                eventTuple().tau_gen_n_gammas.push_back(has_gen_tau ? static_cast<int>(gen_tau.n_gammas) : default_int_value);
                eventTuple().tau_gen_n_gammas_rad.push_back(has_gen_tau ? static_cast<int>(gen_tau.n_gammas_rad)
                                                                : default_int_value);

                eventTuple().tau_decayMode.push_back(tau ? tau->decayMode() : default_int_value);
                eventTuple().tau_oldDecayModeFinding.push_back(tau ? tau->tauID("decayModeFinding") > 0.5f : default_int_value);

                // std::cout<< "reached id fill" << std::endl;
                for(const auto& tau_id_entry : tau_id::GetTauIdDescriptors()) {
                    const auto& desc = tau_id_entry.second;
                    desc.FillTuple(eventTuple, tau, default_value);
                }

                eventTuple().tau_dxy.push_back(tau ? tau->dxy() : default_value);
                eventTuple().tau_dxy_error.push_back(tau ? tau->dxy_error() : default_value);
                eventTuple().tau_ip3d.push_back(tau ? tau->ip3d() : default_value);
                eventTuple().tau_ip3d_error.push_back(tau ? tau->ip3d_error() : default_value);

                const pat::PackedCandidate* leadChargedHadrCand = tau
                        ? dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get())
                        : nullptr;
                eventTuple().tau_dz.push_back(leadChargedHadrCand ? leadChargedHadrCand->dz() : default_value);
                eventTuple().tau_dz_error.push_back(leadChargedHadrCand && leadChargedHadrCand->hasTrackDetails()
                        ? leadChargedHadrCand->dzError() : default_value);

                const auto tauTriggerMatch = triggerDescriptors.matchTriggerObjects(*triggerResults, *triggerObjects,
                    tau_ref_p4, triggerNames.triggerNames(), deltaR2Thr, true, true);
                eventTuple().hlt_accept.push_back(tauTriggerMatch.accept.to_ullong());
                eventTuple().hlt_acceptAndMatch.push_back(tauTriggerMatch.acceptAndMatch.to_ullong());
                
                for(const auto& match_entry : tauTriggerMatch.matchResults) {
                    if(!trig_obj.count(match_entry.second.hltObjIndex)){
                        trig_obj[match_entry.second.hltObjIndex] = match_entry.second;
                    }
                }

                auto l1Tau = MatchL1Taus(tau_ref_p4, *l1Taus, deltaR2Thr, 0);
                eventTuple().l1Tau_pt.push_back(l1Tau ? static_cast<float>(l1Tau->polarP4().pt()) : default_value);
                eventTuple().l1Tau_eta.push_back(l1Tau ? static_cast<float>(l1Tau->polarP4().eta()) : default_value);
                eventTuple().l1Tau_phi.push_back(l1Tau ? static_cast<float>(l1Tau->polarP4().phi()) : default_value);
                eventTuple().l1Tau_mass.push_back(l1Tau ? static_cast<float>(l1Tau->polarP4().mass()) : default_value);
                eventTuple().l1Tau_hwIso.push_back(l1Tau ? l1Tau->hwIso() : default_int_value);
                eventTuple().l1Tau_hwQual.push_back(l1Tau ? l1Tau->hwQual() : default_int_value);
            }

            for(const auto& match_entry : trig_obj) {
                const auto& hlt_obj = triggerObjects->at(match_entry.second.hltObjIndex);
                eventTuple().hltObj_types.push_back(match_entry.second.objType);
                eventTuple().hltObj_pt.push_back(static_cast<float>(hlt_obj.polarP4().pt()));
                eventTuple().hltObj_eta.push_back(static_cast<float>(hlt_obj.polarP4().eta()));
                eventTuple().hltObj_phi.push_back(static_cast<float>(hlt_obj.polarP4().phi()));
                eventTuple().hltObj_mass.push_back(static_cast<float>(hlt_obj.polarP4().mass()));
                eventTuple().hltObj_hasPathName.push_back(match_entry.second.hasPathName.to_ullong());
                eventTuple().hltObj_isBestMatch.push_back(match_entry.second.isBestMatch.to_ullong());
                eventTuple().hltObj_miniAODIndex.push_back(match_entry.second.hltObjIndex);
                const size_t hltObj_index = eventTuple().hltObj_pt.size() - 1;
                for(const std::string& filter : match_entry.second.filters) {
                    eventTuple().filter_hltObj.push_back(hltObj_index);
                    const uint32_t hash = SummaryProducerData::GetData().getFilterHash(filter);
                    eventTuple().filter_hash.push_back(hash);
                }
            }
            // std::cout<< "reached total fill" << std::endl;
            eventTuple.Fill();
            // std::cout << "On to the next" << std::endl;
        }
    }

private:
    const double btagThreshold;
    const bool isMC;
    const std::string triggerProcess;

    edm::EDGetTokenT<GenEventInfoProduct> genEvent_token;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfo_token;
    edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_token;
    edm::EDGetTokenT<pat::TauCollection> taus_token;
    edm::EDGetTokenT<pat::JetCollection> jets_token;
    edm::EDGetTokenT<pat::METCollection> met_token;
    edm::EDGetTokenT<edm::TriggerResults> triggerResults_token;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_token;
    edm::EDGetTokenT<l1t::TauBxCollection> l1Taus_token;

    TriggerDescriptorCollection triggerDescriptors;
    const TupleProducerData& data;
};

} // namespace tau_trigger

#include "FWCore/Framework/interface/MakerMacros.h"
using TauTriggerTupleProducer = tau_trigger::TupleProducer;
DEFINE_FWK_MODULE(TauTriggerTupleProducer);
