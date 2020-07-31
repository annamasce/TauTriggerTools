/*! Produces MvaTuple for tau analysis.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#include <Compression.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TauTriggerTools/Common/interface/GenTruthTools.h"
#include "TauTriggerTools/Studies/interface/MvaTuple.h"

namespace tau_trigger {

struct MvaTupleProducerData {
    using Mutex = MvaTuple::Mutex;
    using LockGuard = std::lock_guard<Mutex>;

    std::unique_ptr<MvaTuple> mvaTuple;

public:
    MvaTupleProducerData(TFile& file)
    {
        mvaTuple = std::make_unique<MvaTuple>("taus", &file, false);
    }
};

class MvaTupleProducer : public edm::stream::EDProducer<edm::GlobalCache<MvaTupleProducerData>> {
public:
    MvaTupleProducer(const edm::ParameterSet& cfg, const MvaTupleProducerData* producerData) :
        isMC(cfg.getParameter<bool>("isMC")),
        genParticles_token(mayConsume<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"))),
        taus_token(consumes<pat::TauCollection>(cfg.getParameter<edm::InputTag>("taus"))),
        data(producerData),
        mvaTuple(*data->mvaTuple)
    {
        produces<bool>();
    }

    static std::unique_ptr<MvaTupleProducerData> initializeGlobalCache(const edm::ParameterSet&)
    {
        TFile& file = edm::Service<TFileService>()->file();
        file.SetCompressionAlgorithm(ROOT::kLZ4);
        file.SetCompressionLevel(4);
        return std::make_unique<MvaTupleProducerData>(file);
    }

    static void globalEndJob(MvaTupleProducerData* data)
    {
        MvaTupleProducerData::LockGuard lock(data->mvaTuple->GetMutex());
        data->mvaTuple->Write();
    }

private:
    static constexpr float default_value = ::tau_trigger::DefaultFillValue<float>();
    static constexpr int default_int_value = ::tau_trigger::DefaultFillValue<int>();

    virtual void produce(edm::Event& event, const edm::EventSetup&) override
    {
        event.put(std::make_unique<bool>(true));

        MvaTupleProducerData::LockGuard lock(data->mvaTuple->GetMutex());

        mvaTuple().run  = event.id().run();
        mvaTuple().lumi = event.id().luminosityBlock();
        mvaTuple().evt  = event.id().event();

        edm::Handle<std::vector<reco::GenParticle>> hGenParticles;
        if(isMC)
            event.getByToken(genParticles_token, hGenParticles);
        auto genParticles = hGenParticles.isValid() ? hGenParticles.product() : nullptr;

        edm::Handle<pat::TauCollection> taus;
        event.getByToken(taus_token, taus);
        for(size_t tau_index = 0; tau_index < taus->size(); ++tau_index) {
            const pat::Tau& tau = taus->at(tau_index);

            mvaTuple().tau_index = tau_index;
            mvaTuple().tau_pt = static_cast<float>(tau.polarP4().pt());
            mvaTuple().tau_eta = static_cast<float>(tau.polarP4().eta());
            mvaTuple().tau_phi = static_cast<float>(tau.polarP4().phi());
            mvaTuple().tau_mass = static_cast<float>(tau.polarP4().mass());
            mvaTuple().tau_charge = tau.charge();

            if(genParticles) {
                const auto gen_match = analysis::gen_truth::LeptonGenMatch(tau.polarP4(), *genParticles);
                mvaTuple().lepton_gen_match = static_cast<int>(gen_match.match);
            } else {
                mvaTuple().lepton_gen_match = default_int_value;
            }

            mvaTuple().tau_dxy_pca_x = tau.dxy_PCA().x();
            mvaTuple().tau_dxy_pca_y = tau.dxy_PCA().y();
            mvaTuple().tau_dxy_pca_z = tau.dxy_PCA().z();
            mvaTuple().tau_dxy = tau.dxy();
            mvaTuple().tau_dxy_error = tau.dxy_error();
            mvaTuple().tau_ip3d = tau.ip3d();
            mvaTuple().tau_ip3d_error = tau.ip3d_error();
            const bool has_sv = tau.hasSecondaryVertex();
            mvaTuple().tau_hasSecondaryVertex = has_sv;
            mvaTuple().tau_sv_x = has_sv ? tau.secondaryVertexPos().x() : default_value;
            mvaTuple().tau_sv_y = has_sv ? tau.secondaryVertexPos().y() : default_value;
            mvaTuple().tau_sv_z = has_sv ? tau.secondaryVertexPos().z() : default_value;
            mvaTuple().tau_flightLength_x = tau.flightLength().x();
            mvaTuple().tau_flightLength_y = tau.flightLength().y();
            mvaTuple().tau_flightLength_z = tau.flightLength().z();
            mvaTuple().tau_flightLength_sig = tau.flightLengthSig();

            mvaTuple.Fill();
        }
    }

private:
    const bool isMC;

    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token;
    edm::EDGetTokenT<pat::TauCollection> taus_token;

    const MvaTupleProducerData* data;
    MvaTuple& mvaTuple;
};

} // namespace tau_trigger

#include "FWCore/Framework/interface/MakerMacros.h"
using TauTriggerMvaTupleProducer = tau_trigger::MvaTupleProducer;
DEFINE_FWK_MODULE(TauTriggerMvaTupleProducer);
