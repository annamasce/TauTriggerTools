#!/usr/bin/env python

import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT

parser = argparse.ArgumentParser(description='Skim full tuple.')
# parser.add_argument('--input', required=True, type=str, nargs='+', help="input files")
parser.add_argument('--config', required=True, type=str, help="config with triggers description")
parser.add_argument('--selection', required=True, type=str, help="tau selection", choices=['DeepTau', 'MVA'])
# parser.add_argument('--output', required=True, type=str, help="output file")
parser.add_argument('--type', required=True, type=str, help="data or mc")
parser.add_argument('--pu', required=False, type=str, default=None,
                    help="file with the pileup profile for the data taking period")
args = parser.parse_args()


#
# Create new branches (as ['var1/F','var2/I',...]) on the given tree
# Returns struct with the new variables
#
import ROOT

cType = {
    'b': 'UChar_t',
    'S': 'Short_t',
    's': 'UShort_t',
    'I': 'Int_t',
    'i': 'UInt_t',
    'F': 'Float_t',
    'D': 'Double_t',
    'L': 'Long64_t',
    'l': 'ULong64_t',
    'O': 'Bool_t',
}
#
# Function to add new branches to a tree
# the ProcessLine makes a new structure newVars that contains a sort of list of all branches as objects, as far as I understand
#
def makeBranches(tree, branches):
    branches = [tuple(branch.split('/')) for branch in branches]        
    ROOT.gROOT.ProcessLine('struct newVars {' + ';'.join([cType[t] + ' ' + name for name, t in branches]) + ';};')
    
    from ROOT import newVars
    newVars = newVars()
    
    for name, t in sorted(branches):
        tree.Branch(name.split('[')[0], ROOT.AddressOf(newVars, name.split('[')[0]), name+ '/' + t)
    return newVars



path_prefix = '' if 'TauTriggerTools' in os.getcwd() else 'TauTriggerTools/'
sys.path.insert(0, path_prefix + 'Common/python')
from TauTriggerTools.Common.AnalysisTypes import *
from TauTriggerTools.Common.AnalysisTools import *
import TauTriggerTools.Common.TriggerConfig as TriggerConfig
ROOT.ROOT.EnableImplicitMT(4)
ROOT.gROOT.SetBatch(True)
# ROOT.gInterpreter.Declare('#include "TauTriggerTools/TauTagAndProbe/interface/PyInterface.h"')

# in_file_path = '/storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_20/src/TauTriggerTools/TauTagAndProbe/test/eventTuple.root'
# in_file_path = '/storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_20/src/TauTriggerTools/TauTagAndProbe/test/eventTuple_WZ.root'
#in_file_path = '/storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_20/src/TauTriggerTools/Studies/data/Factorization/Input/eventTuple_DY.root'
#out_file_path = '/storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_20/src/TauTriggerTools/Studies/data/Factorization/Input/eventTupleSkimmed_DY_'+args.selection+'.root'
#in_file_path = '/storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_20/src/TauTriggerTools/Studies/data/Factorization/Input/eventTuple_WZ.root'
in_file_path = '/storage_mnt/storage/user/lwezenbe/CMSSW_10_2_20/src/TauTriggerTools/TauTagAndProbe/test/eventTuple.root'
out_file_path = '/storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_20/src/TauTriggerTools/Studies/data/Factorization/Input/eventTupleSkimmed_WZ_'+args.selection+'.root'
fIn = ROOT.TFile(in_file_path)
tIn = fIn.Get('events')

fOut = ROOT.TFile(out_file_path, 'RECREATE')

tOut = ROOT.TTree('tOut', 'events')
branches = []
branches.extend(['tau1_pt/F', 'tau1_eta/F', 'tau1_phi/F', 'tau1_mass/F', 'tau1_charge/I', 'tau1_decayMode/I',
    'tau1_byIsolationMVArun2017v2DBoldDMwLT2017/I', 'tau1_byDeepTau2017v2p1VSjet/I'])
branches.extend(['tau2_pt/F', 'tau2_eta/F', 'tau2_phi/F', 'tau2_mass/F', 'tau2_charge/I', 'tau2_decayMode/I',
    'tau2_byIsolationMVArun2017v2DBoldDMwLT2017/I', 'tau2_byDeepTau2017v2p1VSjet/I'])
branches.extend(['passedTrigger/O'])

new_vars = makeBranches(tOut, branches)

    

# if args.type not in ['data', 'mc']:
#     raise RuntimeError("Invalid sample type")

#input_vec = ListToStdVector(in_file_path)
input_vec = in_file_path

#
# Delta phi and R function
#
def deltaPhi(phi1, phi2):
    dphi = phi2-phi1
    if dphi > math.pi:   dphi -= 2.0*math.pi
    if dphi <= -math.pi: dphi += 2.0*math.pi
    return abs(dphi)

def deltaR(eta1, eta2, phi1, phi2):
    return math.sqrt(deltaPhi(phi1, phi2)**2 + (eta1-eta2)**2)

def passFilters(match_desc, hltObj_index, filter_hltObj, filter_hash):
    for filter_ref in match_desc['filter_hashes']:
        filter_found = True
        for n in xrange(len(filter_hltObj)):
            filter_found = filter_hltObj[n] == hltObj_index and filter_hash[n] == filter_ref
            if filter_found: break
        if not filter_found: return False
    
    return True

def tauPassTrigger(chain, index, match_descriptions, channel_index, deltaRThr2, previous_tau_index = None):
    desc_iter = match_descriptions[channel_index]
    for match_desc in desc_iter:
        if (chain.hlt_accept[index] & match_desc['match_mask']) == 0: continue
        if('min_run' in match_desc.keys() and match_desc['min_run'] >= 0 and chain.run < match_desc['min_run']): continue
        if('max_run' in match_desc.keys() and match_desc['max_run'] >= 0 and chain.run >= match_desc['max_run']): continue
        if('l1Tau_pt' in match_desc.keys() and match_desc['l1Tau_pt'] >= 0 and chain.l1Tau_pt < match_desc['l1Tau_pt']): continue
        if('l1Tau_hwIso' in match_desc.keys() and match_desc['l1Tau_hwIso'] >= 0 and chain.l1Tau_hwIso < match_desc['l1Tau_hwIso']): continue

        for n in xrange(len(chain.hltObj_pt[index])):
            if (chain.hltObj_types[index][n] & 4) == 0: continue
            # print chain.hltObj_hasPathName[index][n], match_desc['match_mask']
            if (chain.hltObj_hasPathName[index][n] & match_desc['match_mask']) == 0: continue
            dR2 = deltaR(chain.tau_eta[index], chain.hltObj_eta[index][n], chain.tau_phi[index], chain.hltObj_phi[index][n])
            if dR2 >= deltaRThr2 : continue
            if passFilters(match_desc, n, chain.filter_hltObj[index], chain.filter_hash[index]): return True
        
    return False

# # if args.type == 'mc':
#     # if args.pu is None:
#     #     raise RuntimeError("Pileup file should be provided for mc.")
#     # data_pu_file = ROOT.TFile(args.pu, 'READ')
#     # data_pu = data_pu_file.Get('pileup')
#     # df_all = ROOT.RDataFrame('all_events', input_vec)
#     # mc_pu = df_all.Histo1D(ROOT.RDF.TH1DModel(data_pu), 'npu')
#     # ROOT.PileUpWeightProvider.Initialize(data_pu, mc_pu.GetPtr())

trig_descriptors, channel_triggers = TriggerConfig.Load(args.config)
trigger_dict, filter_dict = TriggerConfig.LoadTriggerDictionarySelf(input_vec)

channels = {}
match_descriptions = {}
for channel_name, channel_trig_descs in channel_triggers.items():
    channel_id = ParseEnum(Channel, channel_name)
    channels[channel_name] = channel_id
    match_descriptions[channel_id] = []
    for desc in channel_trig_descs:
        match_desc = {}
        if 'sample_types' in desc and args.type not in desc['sample_types']: continue
        if desc['leg_types'][-1] != 'tau': continue
        pattern = '^{}.*'.format(desc['name'])
        print pattern
        hlt_paths = TriggerConfig.GetMatchedTriggers(trigger_dict, pattern)
        print hlt_paths
        match_desc['match_mask'] = int(TriggerConfig.GetMatchMask(hlt_paths))
        filter_names = desc['filters'][-1]
        match_desc['filter_hashes'] = ListToStdVector([ filter_dict[f] for f in filter_names ], elem_type='UInt_t')
        if 'min_run' in desc and args.type == 'data':
            matc_desc['min_run'] = desc['min_run']
        if 'max_run' in desc and args.type == 'data':
            max_run['max_run'] = desc['max_run']
        sel_name = 'selection_' + channel_name
        if sel_name in desc:
            if 'hltObj_pt' in desc[sel_name]:
                match_desc['hltObj_pt'] = desc[sel_name]['hltObj_pt']
            if 'l1Tau_pt' in desc[sel_name]:
                match_desc['l1Tau_pt'] = desc[sel_name]['l1Tau_pt']
            if 'l1Tau_hwIso' in desc[sel_name]:
                match_desc['l1Tau_hwIso'] = desc[sel_name]['l1Tau_hwIso']
        match_descriptions[channel_id].append(match_desc)

selection_id = ParseEnum(TauSelection, args.selection)

selection_WP_map = {'VVVLoose' : DiscriminatorWP.VVVLoose,
                    'VVLoose'  : DiscriminatorWP.VVLoose,
                    'VLoose'   : DiscriminatorWP.VLoose,
                    'Loose'    : DiscriminatorWP.Loose,
                    'Medium'   : DiscriminatorWP.Medium,
                    'Tight'    : DiscriminatorWP.Tight,
                    'VTight'    : DiscriminatorWP.VTight,
                    'VVTight'    : DiscriminatorWP.VVTight,
                    'VVVTight'    : DiscriminatorWP.VVVTight,
                    }

from TauTriggerTools.Studies.helpers import progress
deltaRThr2 = 0.5
event_range = xrange(tIn.GetEntries())
for entry in event_range:
    tIn.GetEntry(entry)
    progress(entry, len(event_range))

    if tIn.muon_pt <= 27 or tIn.muon_iso >= 0.1 or tIn.muon_pt <= 30: continue
    tau_indices = []
    #Select 2 taus
    for tau in xrange(len(tIn.tau_pt)):
        if tIn.tau_pt[tau] <= 20 or abs(tIn.tau_eta[tau]) >= 2.1 or tIn.tau_decayMode[tau] == 5 or tIn.tau_decayMode[tau] == 6: continue
        if (tIn.tau_sel[tau] & selection_id) == 0: continue
        if tIn.tau_gen_match == 5: continue
        tau_indices.append(tau)

    if len(tau_indices) != 2: continue

    if tIn.muon_charge + tIn.tau_charge[tau_indices[0]] + tIn.tau_charge[tau_indices[1]] == 0: continue

    both_passed = True
    for tau in tau_indices:
        if not tauPassTrigger(tIn, tau, match_descriptions, 4, deltaRThr2):
            both_passed= False
    

    new_vars.tau1_pt = tIn.tau_pt[tau_indices[0]]
    new_vars.tau1_eta = tIn.tau_eta[tau_indices[0]]
    new_vars.tau1_phi = tIn.tau_phi[tau_indices[0]]
    new_vars.tau1_mass = tIn.tau_mass[tau_indices[0]]
    new_vars.tau1_charge = tIn.tau_charge[tau_indices[0]]
    new_vars.tau1_decayMode = tIn.tau_decayMode[tau_indices[0]]
    new_vars.tau1_byIsolationMVArun2017v2DBoldDMwLT2017 = tIn.byIsolationMVArun2017v2DBoldDMwLT2017[tau_indices[0]]
    new_vars.tau1_byDeepTau2017v2p1VSjet = tIn.byDeepTau2017v2p1VSjet[tau_indices[0]]

    new_vars.tau2_pt = tIn.tau_pt[tau_indices[1]]
    new_vars.tau2_eta = tIn.tau_eta[tau_indices[1]]
    new_vars.tau2_phi = tIn.tau_phi[tau_indices[1]]
    new_vars.tau2_mass = tIn.tau_mass[tau_indices[1]]
    new_vars.tau2_charge = tIn.tau_charge[tau_indices[1]]
    new_vars.tau2_decayMode = tIn.tau_decayMode[tau_indices[1]]
    new_vars.tau2_byIsolationMVArun2017v2DBoldDMwLT2017 = tIn.byIsolationMVArun2017v2DBoldDMwLT2017[tau_indices[1]]
    new_vars.tau2_byDeepTau2017v2p1VSjet = tIn.byDeepTau2017v2p1VSjet[tau_indices[1]]

    new_vars.passedTrigger = both_passed
    
    tOut.Fill()

tOut.AutoSave()
fOut.Close()
  
