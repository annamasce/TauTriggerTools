/*! Definition of a tuple with variables used as inputs to deepTau ID.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#pragma once

#include "TauTriggerTools/Common/interface/SmartTree.h"

#define VAR2(type, name1, name2) VAR(type, name1) VAR(type, name2)
#define VAR3(type, name1, name2, name3) VAR2(type, name1, name2) VAR(type, name3)
#define VAR4(type, name1, name2, name3, name4) VAR3(type, name1, name2, name3) VAR(type, name4)

#define COUNTER_TUPLE_DATA() \
    /* Event Variables */ \
    VAR(UInt_t, run) /* run number */ \
    VAR(UInt_t, lumi) /* lumi section */ \
    VAR(ULong64_t, evt) /* event number */ \
    VAR(Int_t, npv) /* number of primary vertices */ \
    VAR(Float_t, npu) /* number of in-time pu interactions added to the event */ \
    VAR(std::vector<Int_t>, lepton_gen_match) /* gen match info */ \
    VAR(std::vector<Float_t>, gen_tau_pt) /* pt of the gen tau */ \
    VAR(std::vector<Float_t>, gen_tau_eta) /* eta of the gen tau */ \
    VAR(std::vector<Float_t>, gen_tau_phi) /* phi of the gen tau */ \
    VAR(std::vector<Float_t>, gen_tau_e) /* energy of the gen tau */ \
    VAR(std::vector<Float_t>, tau_pt) /* pt of the tau */ \
    VAR(std::vector<Float_t>, tau_eta) /* eta of the tau */ \
    VAR(std::vector<Float_t>, tau_phi) /* phi of the tau */ \
    VAR(std::vector<Float_t>, tau_e) /* energy of the tau */ \
    VAR(std::vector<Float_t>, tau_vz) /* vz of the tau */ \
    VAR(std::vector<Float_t>, deepTau_VSe) /* deepTau_VSe raw discriminator */ \
    VAR(std::vector<Float_t>, deepTau_VSmu) /* deepTau_VSmu raw discriminator */ \
    VAR(std::vector<Float_t>, deepTau_VSjet) /* deepTau_VSjet raw discriminator */ \
    VAR(std::vector<Float_t>, met_pt) /* pt of the MET */ \
    VAR(std::vector<bool>, tau_passedLastFilter) /* bool if the passed the last Filter */ \
    VAR(std::vector<bool>, tau_passedTrackFilter) /* bool if the passed the trackFinding Filter */ \
    VAR(std::vector<Float_t>, gen_met_calo_pt) /* pt of the gen met calo */ \
    VAR(std::vector<Float_t>, gen_met_calo_phi) /* phi of the gen met calo*/ \
    VAR(std::vector<Float_t>, gen_met_true_pt) /* pt of the gen met true */ \
    VAR(std::vector<Float_t>, gen_met_true_phi) /* phi of the gen met true */ \
    VAR(std::vector<Float_t>, l2nn_output) /* output of L2 NN */ \
    VAR(std::vector<Float_t>, l1_pt) /* L1 tau pt */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(counter_tau, TauIdVars, CounterTuple, COUNTER_TUPLE_DATA, "counter_taus")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(counter_tau, CounterTuple, COUNTER_TUPLE_DATA)
#undef VAR
#undef VAR2
#undef VAR3
#undef VAR4
#undef COUNTER_TUPLE_DATA

namespace counter_tau {

template<typename T>
constexpr T DefaultFillValue() { return std::numeric_limits<T>::lowest(); }
template<>
constexpr float DefaultFillValue<float>() { return -999.; }
template<>
constexpr int DefaultFillValue<int>() { return -999; }

} // namespace counter_tau