import os
from helpers_old import makeDirIfNeeded

ptregions = ['high']
ptthreshold = '45'

for ptreg in ptregions:
    makeDirIfNeeded('data/'+ptreg+'/'+ptthreshold)
    #os.system('hadd -f data/'+ptreg+'/'+ptthreshold+'/tauTriggerFactorization2018-2D-v7.root data/Output/'+ptreg+'/'+ptthreshold+'/*root')
    #os.system('rm data/Output/'+ptreg+'/'+ptthreshold+'/*root')
    os.system('hadd -f data/'+ptreg+'/40/tauTriggerFactorization2018-2D-v7.root data/Output/test/'+ptreg+'/40/*root')
    os.system('rm data/Output/test/'+ptreg+'/40/*root')
