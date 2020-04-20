import ROOT
import os
import Sample
import jobSubmitter as sub
from helpers_old import makeDirIfNeeded

sampleList = Sample.createSampleList('/user/lwezenbe/private/CMSSW_10_2_17/src/TauTriggerTools/Factorization/data/inputFiles_2018.conf')
sample = Sample.getSampleFromList(sampleList, 'noTagAndProbe')

ptthreshold = '40'
#ptleading = '45'

for subJob in xrange(sample.splitJobs):
    #for region in ['high', 'low', 'all']:
    #for region in ['high', 'low']:
    for region in ['all']:
        
        log = "/storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_17/src/TauTriggerTools/Factorization/log/log_"+str(subJob)+region+".txt"
        command = 'python /storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_17/src/TauTriggerTools/Factorization/produceFactorization.py --subJob='+str(subJob) + ' --ptRegion=' +region + ' --ptThreshold='+ptthreshold
#        command = 'python /storage_mnt/storage/user/lwezenbe/private/CMSSW_10_2_17/src/TauTriggerTools/Factorization/produceFactorization.py --subJob='+str(subJob) + ' --ptRegion=' +region + ' --ptThreshold='+ptthreshold +' --ptThresholdLeading='+ptleading
        sub.launchCream02(command, log, False, 'Factorization_'+ str(subJob))

