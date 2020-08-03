/*! Various utility functions.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#include "TauTriggerTools/Common/interface/PatHelpers.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TauTriggerTools/Common/interface/AnalysisTypes.h"
#include "TauTriggerTools/Common/interface/CutTools.h"


namespace tau_trigger {

double MuonIsolation(const pat::Muon& muon)
{
    const double pfIso = muon.pfIsolationR04().sumChargedHadronPt
                         + std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt
                                    + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt);
    return pfIso / muon.polarP4().pt();
}

std::vector<TauEntry> CollectTaus(const LorentzVectorM& muon_p4, const pat::TauCollection& taus,
                                  const std::vector<gen_truth::LeptonMatchResult>& genLeptons, double deltaR2Thr)
{
    static const std::string mvaIdName = "byIsolationMVArun2017v2DBoldDMwLTraw2017";
    static const std::string deepIdName = "byDeepTau2017v2p1VSjetraw";
    std::map<TauSelection, const pat::Tau*> best_tau;
    for(const auto& tau : taus) {
        auto leadChargedHadrCand = dynamic_cast<const pat::PackedCandidate*>(tau.leadChargedHadrCand().get());
        if(tau.polarP4().pt() > 18 && std::abs(tau.polarP4().eta()) < 2.3
                && leadChargedHadrCand && std::abs(leadChargedHadrCand->dz()) < 0.2
                && reco::deltaR2(muon_p4, tau.polarP4()) > deltaR2Thr) {
            const bool pass_mva_sel = tau.isTauIDAvailable(mvaIdName) && tau.tauID("againstMuonLoose3") > 0.5f;
            const bool pass_deep_sel = tau.isTauIDAvailable(deepIdName)
                && tau.tauID("byVVVLooseDeepTau2017v2p1VSe") > 0.5f
                && tau.tauID("byVLooseDeepTau2017v2p1VSmu") > 0.5f;
            if((pass_mva_sel || pass_deep_sel) && (!best_tau.count(TauSelection::pt)
                    || best_tau.at(TauSelection::pt)->polarP4().pt() < tau.polarP4().pt()))
                best_tau[TauSelection::pt] = &tau;
            if(pass_mva_sel && (!best_tau.count(TauSelection::MVA)
                    || best_tau.at(TauSelection::MVA)->tauID(mvaIdName) < tau.tauID(mvaIdName)))
                best_tau[TauSelection::MVA] = &tau;
            if(pass_deep_sel && (!best_tau.count(TauSelection::DeepTau)
                    || best_tau.at(TauSelection::DeepTau)->tauID(deepIdName) < tau.tauID(deepIdName)))
                best_tau[TauSelection::DeepTau] = &tau;
        }
    }
    std::map<const pat::Tau*, TauEntry> selected_taus;
    const gen_truth::LeptonMatchResult selected_gen_tau = SelectGenLeg(genLeptons, true);
    const bool has_selected_gen_tau = selected_gen_tau.match != GenLeptonMatch::NoMatch;
    bool selected_gen_tau_stored = false;
    for(const auto& entry : best_tau) {
        const pat::Tau* reco_tau = entry.second;
        if(!selected_taus.count(reco_tau)) {
            const auto gen_tau = gen_truth::LeptonGenMatch(reco_tau->polarP4(), genLeptons);
            const bool has_gen_tau = gen_tau.match != GenLeptonMatch::NoMatch;
            selected_taus[reco_tau] = TauEntry{reco_tau, gen_tau, 0};
            if(has_selected_gen_tau && has_gen_tau
                    && selected_gen_tau.gen_particle_firstCopy == gen_tau.gen_particle_firstCopy) {
                selected_gen_tau_stored = true;
                selected_taus[reco_tau].selection |= static_cast<unsigned>(TauSelection::gen);
            }
        }
        selected_taus[reco_tau].selection |= static_cast<unsigned>(entry.first);
    }
    if(has_selected_gen_tau && !selected_gen_tau_stored) {
        const pat::Tau* reco_tau = nullptr;
        for(const auto& tau : taus) {
            const auto gen_tau = gen_truth::LeptonGenMatch(tau.polarP4(), genLeptons);
            if(gen_tau.match != GenLeptonMatch::NoMatch
                    && gen_tau.gen_particle_firstCopy == selected_gen_tau.gen_particle_firstCopy) {
                reco_tau = &tau;
                break;
            }
        }
        if(selected_taus.count(reco_tau))
            throw exception("Inconsistency in CollectTaus algorithm.");
        selected_taus[reco_tau] = TauEntry{reco_tau, selected_gen_tau, static_cast<unsigned>(TauSelection::gen)};
    }

    std::vector<TauEntry> result;
    for(const auto& entry : selected_taus)
        result.push_back(entry.second);
    return result;
}


bool IsGoodBaselineTau(const pat::Tau& tau, double deltaR2Thr){
    auto leadChargedHadrCand = dynamic_cast<const pat::PackedCandidate*>(tau.leadChargedHadrCand().get());
    return (tau.polarP4().pt() > 18 && std::abs(tau.polarP4().eta()) < 2.3
            && leadChargedHadrCand && std::abs(leadChargedHadrCand->dz()) < 0.2);
}

bool IsBetterTauPair(std::vector<const pat::Tau*>& tau_pair_1, const std::vector<const pat::Tau*>& tau_pair_2, const std::string IdName){
    if(tau_pair_1[0]->tauID(IdName) < tau_pair_2[0]->tauID(IdName)){
        return true;
    } else if(tau_pair_1[0]->tauID(IdName) == tau_pair_2[0]->tauID(IdName)) {
        if(tau_pair_1[0]->pt() < tau_pair_2[0]->pt()) {
            return true;
        } else if(tau_pair_1[0]->pt() == tau_pair_2[0]->pt()) {
            if(tau_pair_1[1]->tauID(IdName) < tau_pair_2[1]->tauID(IdName)){
                return true;
            } else if(tau_pair_1[1]->tauID(IdName) == tau_pair_2[1]->tauID(IdName)){
                if(tau_pair_1[1]->pt() < tau_pair_2[1]->pt()) {
                    return true;
                }else{
                    return false;
                }
            } else {
                return false;
            }
        }else{
            return false;
        }
    }else{
        return false;
        
    }

}

std::vector<std::vector<TauEntry>> CollectTauPairs(const pat::TauCollection& taus,
                                  const std::vector<gen_truth::LeptonMatchResult>& genLeptons, double deltaR2Thr, cuts::Cutter<>& cut)
{
    static const std::string mvaIdName = "byIsolationMVArun2017v2DBoldDMwLTraw2017";
    static const std::string deepIdName = "byDeepTau2017v2p1VSjetraw";
    std::map<TauSelection, std::vector<const pat::Tau*>> best_tau_pair;
    int tau1_index = 0;
    bool has_any_pair = taus.size() > 1;
    bool has_pair = false;
    for(const auto& tau1 : taus) {
        tau1_index++;
        if(!IsGoodBaselineTau(tau1, deltaR2Thr)) continue;

        const bool pass_mva_sel_tau1 = tau1.isTauIDAvailable(mvaIdName) && tau1.tauID("againstMuonLoose3") > 0.5f;
        const bool pass_deep_sel_tau1 = tau1.isTauIDAvailable(deepIdName)
            && tau1.tauID("byVVVLooseDeepTau2017v2p1VSe") > 0.5f
            && tau1.tauID("byVLooseDeepTau2017v2p1VSmu") > 0.5f;
        int tau2_index = 0;
        for(const auto& tau2 : taus) {
            tau2_index++;
            if (tau2_index <= tau1_index)       continue;
            if(!IsGoodBaselineTau(tau2, deltaR2Thr)) continue;

            const bool pass_mva_sel_tau2 = tau2.isTauIDAvailable(mvaIdName) && tau2.tauID("againstMuonLoose3") > 0.5f;
            const bool pass_deep_sel_tau2 = tau2.isTauIDAvailable("byDeepTau2017v2p1VSjetraw")
                && tau2.tauID("byVVVLooseDeepTau2017v2p1VSe") > 0.5f
                && tau2.tauID("byVLooseDeepTau2017v2p1VSmu") > 0.5f;

            const float dR = deltaR (tau1, tau2);
            if(dR < 0.5) continue; 

            if(pass_deep_sel_tau1 && pass_deep_sel_tau2){
                std::vector<const pat::Tau*> tau_pair;
                if(tau1.tauID(deepIdName) > tau2.tauID(deepIdName)){
                    tau_pair.push_back(&tau1);
                    tau_pair.push_back(&tau2);
                }else{
                    tau_pair.push_back(&tau2);
                    tau_pair.push_back(&tau1);
                }
                if(!best_tau_pair.count(TauSelection::DeepTau)){
                    has_pair = true;
                    best_tau_pair[TauSelection::DeepTau] = tau_pair;
                } else if(IsBetterTauPair(best_tau_pair[TauSelection::DeepTau], tau_pair, deepIdName)) {
                    best_tau_pair[TauSelection::DeepTau] = tau_pair;
                }

            }

            if(pass_mva_sel_tau1 && pass_mva_sel_tau2){
                std::vector<const pat::Tau*> tau_pair;
                if(tau1.tauID(mvaIdName) > tau2.tauID(mvaIdName)){
                    tau_pair.push_back(&tau1);
                    tau_pair.push_back(&tau2);
                }else{
                    tau_pair.push_back(&tau2);
                    tau_pair.push_back(&tau1);
                }
                if(!best_tau_pair.count(TauSelection::MVA)){
                    has_pair = true;
                    best_tau_pair[TauSelection::MVA] = tau_pair;
                } else if(IsBetterTauPair(best_tau_pair[TauSelection::MVA], tau_pair, mvaIdName)) {
                    best_tau_pair[TauSelection::MVA] = tau_pair;
                }
            }
        }
    }

    cut(has_any_pair ,"has_any_tau_pair");
    cut(has_pair, "has_selected_tau_pair");


    // Translate taus into TauEntries
    std::map<const pat::Tau*, TauEntry> selected_taus;
    for(const auto& entry : best_tau_pair){
        const pat::Tau* reco_tau1 = entry.second[0];
        const pat::Tau* reco_tau2 = entry.second[1];
        if(!selected_taus.count(reco_tau1)) {
            const auto gen_tau = gen_truth::LeptonGenMatch(reco_tau1->polarP4(), genLeptons);
            selected_taus[reco_tau1] = TauEntry{reco_tau1, gen_tau, 0};
        }
        selected_taus[reco_tau1].selection |= static_cast<unsigned>(entry.first);
 
        if(!selected_taus.count(reco_tau2)) {
            const auto gen_tau = gen_truth::LeptonGenMatch(reco_tau2->polarP4(), genLeptons);
            selected_taus[reco_tau2] = TauEntry{reco_tau2, gen_tau, 0};
        }
        selected_taus[reco_tau2].selection |= static_cast<unsigned>(entry.first);
    }

    // Gather pairs of TauEntries again
    std::map<std::vector<const pat::Tau*>, std::vector<TauEntry>> selected_tau_pairs;

    for(const auto& entry : best_tau_pair){

        if(selected_tau_pairs.empty()){
            std::vector<TauEntry> pair;
            pair.push_back(selected_taus[entry.second[0]]);
            pair.push_back(selected_taus[entry.second[1]]);
            selected_tau_pairs[entry.second] = pair; 
        } else {
            bool already_contained = false;
            // for(selected_entry = selected_tau_pairs.begin(); selected_entry != selected_tau_pairs.end(); selected_entry++){
            for(auto const& selected_entry : selected_tau_pairs){
                if(std::count(selected_entry.first.begin(), selected_entry.first.end(), entry.second[0]) and std::count(selected_entry.first.begin(), selected_entry.first.end(), entry.second[1])){
                    already_contained = true;
                }
            }
            if(!already_contained){
                std::vector<TauEntry> pair;
                pair.push_back(selected_taus[entry.second[0]]);
                pair.push_back(selected_taus[entry.second[1]]);
                selected_tau_pairs[entry.second] = pair; 
            }
        }
    }    
    std::vector<std::vector<TauEntry>> result;
    for(const auto& entry : selected_tau_pairs){
        result.push_back(entry.second);
    }
    return result;
}

bool PassBtagVeto(const LorentzVectorM& tau_p4,
                  const pat::JetCollection& jets, double btagThreshold, double deltaR2Thr)
{
    if(btagThreshold > 0) {
        for(const pat::Jet& jet : jets) {
            const auto btag = jet.bDiscriminator("pfDeepFlavourJetTags:probb")
                              + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")
                              + jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
            if(reco::deltaR2(tau_p4, jet.polarP4()) > deltaR2Thr
                    && jet.polarP4().pt() > 20 && std::abs(jet.polarP4().eta()) < 2.4
                    && btag > btagThreshold)
                return false;
        }
    }
    return true;
}

gen_truth::LeptonMatchResult SelectGenLeg(const std::vector<gen_truth::LeptonMatchResult>& genLeptons, bool is_tau)
{
    static const std::map<bool, std::set<GenLeptonMatch>> all_matches = {
        { true, { GenLeptonMatch::Tau } },
        { false, { GenLeptonMatch::Muon, GenLeptonMatch::TauMuon } },
    };
    const auto& matches = all_matches.at(is_tau);
    gen_truth::LeptonMatchResult leg;
    for(const auto& lepton : genLeptons) {
        if(matches.count(lepton.match) && (leg.match == GenLeptonMatch::NoMatch
                    || leg.visible_p4.pt() < lepton.visible_p4.pt())) {
            leg = lepton;
        }
    }
    return leg;
}

} // namespace tau_trigger
