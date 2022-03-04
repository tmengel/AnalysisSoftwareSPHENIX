input=/home/tmengel/sPHENIX/TannerEval/macro/treeAnalysis_sphenix/G4sPHENIX_TannerEval.root
geometry=/home/tmengel/sPHENIX/TannerEval/macro/treeAnalysis_sphenix/geometry.root



root -x -l -b -q 'treeProcessing_simple_sphenix.C("'$input'","'$geometry'","test_out")'
