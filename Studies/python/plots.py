import os

import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--selection',             action='store',       required=True, choices=['DeepTau', 'MVA']) 

args = argParser.parse_args()
ptregions = ['high']
#var = ['2D', 'leadingpt', 'subleadingpt']
var = ['leadingpt', 'subleadingpt']
#var = ['2D']
wp= ['vloose', 'loose', 'medium', 'tight', 'vtight']
ptthreshold = '40'
for ptreg in ptregions:
    for v in var:
        for i in wp:
	    print 'python producePlots.py --ptRegion='+ ptreg +' --WP='+i+' --var='+v+ ' --ptThreshold='+ptthreshold + ' --selection='+args.selection
            os.system('python producePlots.py --ptRegion='+ ptreg +' --WP='+i+' --var='+v+ ' --ptThreshold='+ptthreshold + ' --selection='+args.selection)
