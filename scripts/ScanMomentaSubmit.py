import os
import sys
import numpy as np
#import ParameterVariations
## Using this version to only do the extra variations
import ParameterVariationsExtra


## Flag to submit jobs or not
SubmitOrNot = False

## Folder for storage (this path)
GENERAL_SCRATCH = os.path.dirname(os.path.realpath(__file__))

## READ IN MOMENTA FOR NUCLEI PER PID
nuclei = ["c","o","al","pb","fe","cu"]
pions = {'piP' : 211, 'piM': -211}

for nucleus in nuclei:    
    for pion in pions:

        ## Merge into single file for easiness
        os.system("cat ../pion_data/%s_*_%s*.csv  > /tmp/all_%s_%s.csv" % (nucleus,pion,nucleus,pion))

        ## Read in first row (momenta values) into a list
        ## and get the unique elements (don't want to repeat)
        momenta = np.unique(np.genfromtxt("/tmp/all_%s_%s.csv" % (nucleus,pion), delimiter=' ')[:,0])
        
        print "%s momenta %s: " % (nucleus,pion)
        print momenta

        ## For each momenta launch the subroutine to create jobs
        for mom in momenta:
            ParameterVariationsExtra.build_pbs_files("single",
                                                              nucleus,
                                                              mom,
                                                              pions[pion],
                                                              601,
                                                              SubmitOrNot)
