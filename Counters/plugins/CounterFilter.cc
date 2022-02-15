/*! Apply tau trigger selection vetoes.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "TauTriggerTools/Common/interface/AnalysisTypes.h"
#include "TauTriggerTools/Common/interface/CutTools.h"
#include "TauTriggerTools/Common/interface/PatHelpers.h"
#include "TauTriggerTools/Common/interface/GenTruthTools.h"
#include "TauTriggerTools/Common/interface/GenLepton.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "RecoTauTag/RecoTau/interface/DeepTauBase.h"
#include "TauTriggerTools/Counters/interface/CounterTuple.h"

class CounterFilter : public edm::EDFilter {
public:
    // using TauDiscriminator = reco::PFTauDiscriminator;
    using TauDiscriminatorContainer = reco::TauDiscriminatorContainer;

    CounterFilter(const edm::ParameterSet& cfg) :
        isMC(cfg.getParameter<bool>("isMC")),
        gen_level(cfg.getParameter<bool>("gen_level")),
        use_deepTau(cfg.getParameter<bool>("use_deepTau")),
        use_L2NN(cfg.getParameter<bool>("use_L2NN")),
        store_MET(cfg.getParameter<bool>("store_MET")),
        position(cfg.getParameter<std::string>("position")),
        deepTauVSe_inputToken(mayConsume<TauDiscriminatorContainer>(cfg.getParameter<edm::InputTag>("deepTauVSe"))),
        deepTauVSmu_inputToken(mayConsume<TauDiscriminatorContainer>(cfg.getParameter<edm::InputTag>("deepTauVSmu"))),
        deepTauVSjet_inputToken(mayConsume<TauDiscriminatorContainer>(cfg.getParameter<edm::InputTag>("deepTauVSjet"))),
        decayMode_token(consumes<reco::PFTauDiscriminator>(cfg.getParameter<edm::InputTag>("decayModeFindingNewDM"))),
        L2NNoutput_token(mayConsume<std::vector<float>>(cfg.getParameter<edm::InputTag>("L2NNoutput"))),
        l1taus_token(mayConsume<trigger::TriggerFilterObjectWithRefs>(cfg.getParameter<edm::InputTag>("l1taus"))),
        original_taus_token(mayConsume<std::vector<reco::PFTau>>(cfg.getParameter<edm::InputTag>("original_taus"))),
        taus_token(mayConsume<trigger::TriggerFilterObjectWithRefs>(cfg.getParameter<edm::InputTag>("taus"))),
        track_taus_token(mayConsume<trigger::TriggerFilterObjectWithRefs>(cfg.getParameter<edm::InputTag>("track_taus"))),
        MET_token(mayConsume<std::vector<reco::CaloMET>>(cfg.getParameter<edm::InputTag>("MET"))),
        puInfo_token(mayConsume<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("puInfo"))),
        vertices_token(mayConsume<std::vector<reco::Vertex> >(cfg.getParameter<edm::InputTag>("vertices"))),
        genParticles_token(mayConsume<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"))),
        genMETcalo_token(mayConsume<std::vector<reco::GenMET>>(cfg.getParameter<edm::InputTag>("genMETcalo"))),
        genMETtrue_token(mayConsume<std::vector<reco::GenMET>>(cfg.getParameter<edm::InputTag>("genMETtrue")))
    {
        std::string full_name = position+"_counter";
        counterTuple = std::make_shared<counter_tau::CounterTuple>(full_name, &edm::Service<TFileService>()->file(), false);
    }

private:
    static constexpr int default_int_value = ::counter_tau::DefaultFillValue<int>();
    static constexpr float default_value = ::counter_tau::DefaultFillValue<float>();

    virtual bool filter(edm::Event& event, const edm::EventSetup&) override
    {
        bool result = true;

        if(gen_level){
            (*counterTuple)().run  = event.id().run();
            (*counterTuple)().lumi = event.id().luminosityBlock();
            (*counterTuple)().evt  = event.id().event();

            if(isMC){
                edm::Handle<std::vector<reco::GenParticle>> hGenParticles;
                event.getByToken(genParticles_token, hGenParticles);
                genParticles = hGenParticles.isValid() ? hGenParticles.product() : nullptr;

                std::vector<reco_tau::gen_truth::GenLepton> lepton_results = reco_tau::gen_truth::GenLepton::fromGenParticleCollection(*genParticles);

                for(unsigned n = 0; n < lepton_results.size(); ++n){
                    const auto gen_match = lepton_results.at(n);
                    (*counterTuple)().lepton_gen_match.push_back(static_cast<int>(gen_match.kind()));
                    (*counterTuple)().gen_tau_charge.push_back(static_cast<int>(gen_match.charge()));
                    (*counterTuple)().gen_tau_pt.push_back(static_cast<float>(gen_match.visibleP4().pt()));
                    (*counterTuple)().gen_tau_eta.push_back(static_cast<float>(gen_match.visibleP4().eta()));
                    (*counterTuple)().gen_tau_phi.push_back(static_cast<float>(gen_match.visibleP4().phi()));
                    (*counterTuple)().gen_tau_e.push_back(static_cast<float>(gen_match.visibleP4().e()));
                    (*counterTuple)().gen_tau_pt_rad.push_back(static_cast<float>(gen_match.radiatedP4().pt()));
                    (*counterTuple)().gen_tau_eta_rad.push_back(static_cast<float>(gen_match.radiatedP4().eta()));
                    (*counterTuple)().gen_tau_phi_rad.push_back(static_cast<float>(gen_match.radiatedP4().phi()));
                    (*counterTuple)().gen_tau_e_rad.push_back(static_cast<float>(gen_match.radiatedP4().e()));
                    (*counterTuple)().gen_tau_nChargedHadrons.push_back(static_cast<int>(gen_match.nChargedHadrons()));
                    (*counterTuple)().gen_tau_nNeutralHadrons.push_back(static_cast<int>(gen_match.nNeutralHadrons()));
                    (*counterTuple)().gen_tau_nFinalStateElectrons.push_back(static_cast<int>(gen_match.nFinalStateElectrons()));
                    (*counterTuple)().gen_tau_nFinalStateMuons.push_back(static_cast<int>(gen_match.nFinalStateMuons()));

                }

                edm::Handle<std::vector<reco::GenMET>> hGenMETcalo;
                event.getByToken(genMETcalo_token, hGenMETcalo);
                genMETcalo = hGenMETcalo.isValid() ? hGenMETcalo.product() : nullptr;
                if (genMETcalo){
                    (*counterTuple)().gen_met_calo_pt.push_back(static_cast<float>((*genMETcalo).at(0).pt()));
                    (*counterTuple)().gen_met_calo_phi.push_back(static_cast<float>((*genMETcalo).at(0).phi()));
                }

                edm::Handle<std::vector<reco::GenMET>> hGenMETtrue;
                event.getByToken(genMETtrue_token, hGenMETtrue);
                genMETtrue = hGenMETtrue.isValid() ? hGenMETtrue.product() : nullptr;
                if (genMETtrue){
                    (*counterTuple)().gen_met_true_pt.push_back(static_cast<float>((*genMETtrue).at(0).pt()));
                    (*counterTuple)().gen_met_true_phi.push_back(static_cast<float>((*genMETtrue).at(0).phi()));
                }
            }

            counterTuple->Fill();
        }
        else{

            edm::Handle<std::vector<reco::Vertex>> vertices;
            event.getByToken(vertices_token, vertices);
            (*counterTuple)().npv = static_cast<int>(vertices->size());

            edm::Handle<TauDiscriminatorContainer> deepTau_VSe;
            edm::Handle<TauDiscriminatorContainer> deepTau_VSmu;
            edm::Handle<TauDiscriminatorContainer> deepTau_VSjet;
            if(use_deepTau){
                event.getByToken(deepTauVSe_inputToken, deepTau_VSe);
                event.getByToken(deepTauVSmu_inputToken, deepTau_VSmu);
                event.getByToken(deepTauVSjet_inputToken, deepTau_VSjet);
            }

            edm::Handle<reco::PFTauDiscriminator> decayModesNew;
            event.getByToken(decayMode_token, decayModesNew);

            edm::Handle<std::vector<float>> L2NNoutput;
            l1t::TauVectorRef l1Taus;
            auto const& l1TriggeredTaus = event.get(l1taus_token);
            l1TriggeredTaus.getObjects(trigger::TriggerL1Tau, l1Taus);
            if(use_L2NN) event.getByToken(L2NNoutput_token, L2NNoutput);
            for(size_t i=0; i<l1Taus.size(); i++){
                if (use_L2NN) {
                    (*counterTuple)().l2nn_output.push_back(static_cast<float>(L2NNoutput->at(i)));
                }
                (*counterTuple)().l1_pt.push_back(static_cast<float>(l1Taus[i]->pt()));
                (*counterTuple)().l1_eta.push_back(static_cast<float>(l1Taus[i]->eta()));
                (*counterTuple)().l1_phi.push_back(static_cast<float>(l1Taus[i]->phi()));
                (*counterTuple)().l1_hwIso.push_back(static_cast<float>(l1Taus[i]->hwIso()));
            }
            
            edm::Handle<std::vector<reco::CaloMET>> met;
            if(store_MET){
                event.getByToken(MET_token, met);
                for(size_t met_index = 0; met_index < met->size(); ++met_index) {
                    const reco::CaloMET& met_cand = met->at(met_index);
                    (*counterTuple)().met_pt.push_back(static_cast<float>(met_cand.polarP4().pt()));
                }

            }
            else{
                (*counterTuple)().met_pt.push_back(default_value);
            }

            edm::Handle<std::vector<reco::PFTau>> original_taus;
            event.getByToken(original_taus_token, original_taus);

            edm::Handle<trigger::TriggerFilterObjectWithRefs> finalTaus;
            event.getByToken(taus_token, finalTaus);

            trigger::VRpftau tauCandRefVec;
            finalTaus->getObjects(trigger::TriggerTau, tauCandRefVec);

            edm::Handle<trigger::TriggerFilterObjectWithRefs> trackTaus;
            event.getByToken(track_taus_token, trackTaus);

            trigger::VRpftau tauCandRefVec_track;
            finalTaus->getObjects(trigger::TriggerTau, tauCandRefVec_track);

            edm::Handle<std::vector<reco::GenParticle>> hGenParticles;
            if(isMC) {
                edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
                event.getByToken(puInfo_token, puInfo);
                (*counterTuple)().npu = analysis::gen_truth::GetNumberOfPileUpInteractions(puInfo);

                event.getByToken(genParticles_token, hGenParticles);

            }

            genParticles = hGenParticles.isValid() ? hGenParticles.product() : nullptr;

            (*counterTuple)().run  = event.id().run();
            (*counterTuple)().lumi = event.id().luminosityBlock();
            (*counterTuple)().evt  = event.id().event();

            edm::Handle<std::vector<reco::GenMET>> hGenMETcalo;
                event.getByToken(genMETcalo_token, hGenMETcalo);
                genMETcalo = hGenMETcalo.isValid() ? hGenMETcalo.product() : nullptr;

            if (genMETcalo){
                (*counterTuple)().gen_met_calo_pt.push_back(static_cast<float>((*genMETcalo).at(0).pt()));
                (*counterTuple)().gen_met_calo_phi.push_back(static_cast<float>((*genMETcalo).at(0).phi()));
            }
            else{
                (*counterTuple)().gen_met_calo_pt.push_back(default_value);
                (*counterTuple)().gen_met_calo_phi.push_back(default_value);
            }

            edm::Handle<std::vector<reco::GenMET>> hGenMETtrue;
                event.getByToken(genMETtrue_token, hGenMETtrue);
                genMETtrue = hGenMETtrue.isValid() ? hGenMETtrue.product() : nullptr;

            if (genMETtrue){
                (*counterTuple)().gen_met_true_pt.push_back(static_cast<float>((*genMETtrue).at(0).pt()));
                (*counterTuple)().gen_met_true_phi.push_back(static_cast<float>((*genMETtrue).at(0).phi()));
            }
            else{
                (*counterTuple)().gen_met_true_pt.push_back(default_value);
                (*counterTuple)().gen_met_true_phi.push_back(default_value);
            }

            std::vector<reco_tau::gen_truth::GenLepton> lepton_results;            
            if(genParticles){
                lepton_results = reco_tau::gen_truth::GenLepton::fromGenParticleCollection(*genParticles);
            }

            for(size_t orig_tau_index = 0; orig_tau_index < original_taus->size(); ++orig_tau_index) {
                const reco::PFTau& original_tau = original_taus->at(orig_tau_index);
                edm::Ref<reco::PFTauCollection> tauRef(original_taus, orig_tau_index);

                if(genParticles){
                    bool matched = false;
                    int match_id = 6;
                    reco_tau::gen_truth::GenLepton gen_match;
                    double deltaR_thr = 0.2;
                    for(size_t gen_index = 0; gen_index < lepton_results.size(); ++gen_index){
                        const auto total_vis_p4 = lepton_results.at(gen_index).visibleP4() + lepton_results.at(gen_index).radiatedP4();
                        const double deltaR = ROOT::Math::VectorUtil::DeltaR(original_tau.polarP4(),total_vis_p4);
                        if(deltaR < deltaR_thr){
                            matched = true;
                            gen_match = lepton_results.at(gen_index);
                            deltaR_thr = deltaR;
                        }
                    } // loop over gen leptons
                    if(matched){
                        (*counterTuple)().lepton_gen_match.push_back(static_cast<int>(gen_match.kind()));
                        (*counterTuple)().gen_tau_charge.push_back(static_cast<int>(gen_match.charge()));
                        (*counterTuple)().gen_tau_pt.push_back(static_cast<float>(gen_match.visibleP4().pt()));
                        (*counterTuple)().gen_tau_eta.push_back(static_cast<float>(gen_match.visibleP4().eta()));
                        (*counterTuple)().gen_tau_phi.push_back(static_cast<float>(gen_match.visibleP4().phi()));
                        (*counterTuple)().gen_tau_e.push_back(static_cast<float>(gen_match.visibleP4().e()));
                        (*counterTuple)().gen_tau_pt_rad.push_back(static_cast<float>(gen_match.radiatedP4().pt()));
                        (*counterTuple)().gen_tau_eta_rad.push_back(static_cast<float>(gen_match.radiatedP4().eta()));
                        (*counterTuple)().gen_tau_phi_rad.push_back(static_cast<float>(gen_match.radiatedP4().phi()));
                        (*counterTuple)().gen_tau_e_rad.push_back(static_cast<float>(gen_match.radiatedP4().e()));
                        (*counterTuple)().gen_tau_nChargedHadrons.push_back(static_cast<int>(gen_match.nChargedHadrons()));
                        (*counterTuple)().gen_tau_nNeutralHadrons.push_back(static_cast<int>(gen_match.nNeutralHadrons()));
                        (*counterTuple)().gen_tau_nFinalStateElectrons.push_back(static_cast<int>(gen_match.nFinalStateElectrons()));
                        (*counterTuple)().gen_tau_nFinalStateMuons.push_back(static_cast<int>(gen_match.nFinalStateMuons()));
                    }
                    else{
                        (*counterTuple)().lepton_gen_match.push_back(match_id);
                        (*counterTuple)().gen_tau_charge.push_back(0.);
                        (*counterTuple)().gen_tau_pt.push_back(0.);
                        (*counterTuple)().gen_tau_eta.push_back(0.);
                        (*counterTuple)().gen_tau_phi.push_back(0.);
                        (*counterTuple)().gen_tau_e.push_back(0.);
                        (*counterTuple)().gen_tau_pt_rad.push_back(0.);
                        (*counterTuple)().gen_tau_eta_rad.push_back(0.);
                        (*counterTuple)().gen_tau_phi_rad.push_back(0.);
                        (*counterTuple)().gen_tau_e_rad.push_back(0.);
                        (*counterTuple)().gen_tau_nChargedHadrons.push_back(0);
                        (*counterTuple)().gen_tau_nNeutralHadrons.push_back(0);
                        (*counterTuple)().gen_tau_nFinalStateElectrons.push_back(0);
                        (*counterTuple)().gen_tau_nFinalStateMuons.push_back(0);
                    }                    
                } else {
                    (*counterTuple)().lepton_gen_match.push_back(default_int_value);
                    (*counterTuple)().gen_tau_charge.push_back(default_int_value);
                    (*counterTuple)().gen_tau_pt.push_back(default_value);
                    (*counterTuple)().gen_tau_eta.push_back(default_value);
                    (*counterTuple)().gen_tau_phi.push_back(default_value);
                    (*counterTuple)().gen_tau_e.push_back(default_value);
                    (*counterTuple)().gen_tau_pt_rad.push_back(default_value);
                    (*counterTuple)().gen_tau_eta_rad.push_back(default_value);
                    (*counterTuple)().gen_tau_phi_rad.push_back(default_value);
                    (*counterTuple)().gen_tau_e_rad.push_back(default_value);
                    (*counterTuple)().gen_tau_nChargedHadrons.push_back(default_int_value);
                    (*counterTuple)().gen_tau_nNeutralHadrons.push_back(default_int_value);
                    (*counterTuple)().gen_tau_nFinalStateElectrons.push_back(default_int_value);
                    (*counterTuple)().gen_tau_nFinalStateMuons.push_back(default_int_value);
                }

                (*counterTuple)().tau_pt.push_back(static_cast<float>(original_tau.polarP4().pt()));
                (*counterTuple)().tau_eta.push_back(static_cast<float>(original_tau.polarP4().eta()));
                (*counterTuple)().tau_phi.push_back(static_cast<float>(original_tau.polarP4().phi()));
                (*counterTuple)().tau_e.push_back(static_cast<float>(original_tau.polarP4().e()));
                (*counterTuple)().tau_vz.push_back(static_cast<float>(original_tau.vz()));
                (*counterTuple)().tau_decayModeFindingNewDMs.push_back(decayModesNew->value(orig_tau_index));

                if(use_deepTau){
                    (*counterTuple)().deepTau_VSe.push_back(static_cast<float>((*deepTau_VSe)[tauRef].rawValues.at(0)));
                    (*counterTuple)().deepTau_VSmu.push_back(static_cast<float>((*deepTau_VSmu)[tauRef].rawValues.at(0)));
                    (*counterTuple)().deepTau_VSjet.push_back(static_cast<float>((*deepTau_VSjet)[tauRef].rawValues.at(0)));
                }
                else{
                    (*counterTuple)().deepTau_VSe.push_back(default_value);
                    (*counterTuple)().deepTau_VSmu.push_back(default_value);
                    (*counterTuple)().deepTau_VSjet.push_back(default_value);
                }

                bool passed_lastFilter = false;
                for (unsigned int iTau = 0; iTau < tauCandRefVec.size(); iTau++) {
                    const double deltaR = ROOT::Math::VectorUtil::DeltaR(original_tau.polarP4(),tauCandRefVec[iTau]->p4());
                    if(deltaR < 0.01){
                        passed_lastFilter = true;
                        break;
                    }
                }

                (*counterTuple)().tau_passedLastFilter.push_back(passed_lastFilter);

                bool passed_trackFilter = false;
                for (unsigned int iTau = 0; iTau < tauCandRefVec_track.size(); iTau++) {
                    const double deltaR = ROOT::Math::VectorUtil::DeltaR(original_tau.polarP4(),tauCandRefVec_track[iTau]->p4());
                    if(deltaR < 0.01){
                        passed_trackFilter = true;
                        break;
                    }
                }

                (*counterTuple)().tau_passedTrackFilter.push_back(passed_trackFilter);

            } // loop over original taus

            counterTuple->Fill();
        }

        return result;
    }

    void endJob()
    {
        counterTuple->Write();
    }

private:
    const bool isMC, gen_level, use_deepTau, use_L2NN, store_MET;
    std::string position;
    const edm::EDGetTokenT<TauDiscriminatorContainer> deepTauVSe_inputToken;
    const edm::EDGetTokenT<TauDiscriminatorContainer> deepTauVSmu_inputToken;
    const edm::EDGetTokenT<TauDiscriminatorContainer> deepTauVSjet_inputToken;
    const edm::EDGetTokenT<reco::PFTauDiscriminator> decayMode_token;
    const edm::EDGetTokenT<std::vector<float>> L2NNoutput_token;
    edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> l1taus_token;
    edm::EDGetTokenT<std::vector<reco::PFTau>> original_taus_token;
    edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> taus_token;
    edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> track_taus_token;
    edm::EDGetTokenT<std::vector<reco::CaloMET>> MET_token;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfo_token;
    edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_token;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token;
    edm::EDGetTokenT<std::vector<reco::GenMET>> genMETcalo_token;
    edm::EDGetTokenT<std::vector<reco::GenMET>> genMETtrue_token;
    const std::vector<reco::GenMET>* genMETcalo;
    const std::vector<reco::GenMET>* genMETtrue;
    const std::vector<reco::GenParticle>* genParticles;
    std::shared_ptr<TH1F> counter;
    std::shared_ptr<counter_tau::CounterTuple> counterTuple;

};


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CounterFilter);