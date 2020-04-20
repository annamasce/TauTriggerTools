import ROOT, socket, os, shutil, subprocess, time
from math import pi, sqrt
from operator import mul
from numpy import loadtxt
import sys

#
# Check if valid ROOT file exists
#
def isValidRootFile(fname):
  if not os.path.exists(os.path.expandvars(fname)): return False
  if 'pnfs' in fname: fname = 'root://maite.iihe.ac.be'+ fname         #faster for pnfs files + avoids certain unstable problems I had with input/output errors
  f = ROOT.TFile.Open(fname)
  if not f: return False
  try:
    return not (f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered) or f.GetListOfKeys().IsEmpty())
  finally:
    f.Close()

#
# Get object (e.g. hist) from file using key, and keep in memory after closing
#
def getObjFromFile(fname, hname):
    assert isValidRootFile(fname)

    if 'pnfs' in fname: fname = 'root://maite.iihe.ac.be'+ fname         #faster for pnfs file
    try:
      f = ROOT.TFile.Open(fname)
      f.cd()
      htmp = f.Get(hname)
      if not htmp: return None
      ROOT.gDirectory.cd('PyROOT:/')
      res = htmp.Clone()
      return res
    finally:
      f.Close()
