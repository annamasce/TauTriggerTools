import os

ptregions = ['high']
#var = ['2D', 'leadingpt', 'subleadingpt']
var = ['leadingpt', 'subleadingpt']
#var = ['2D']
wp= ['vlooseTauMVA', 'looseTauMVA', 'mediumTauMVA', 'tightTauMVA', 'vtightTauMVA']
ptthreshold = '40'
for ptreg in ptregions:
    for v in var:
        for i in wp:
            os.system('python producePlots.py --ptRegion='+ ptreg +' --WP='+i+' --var='+v+ ' --ptThreshold='+ptthreshold)
