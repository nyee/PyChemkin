import logging
logging.disable(logging.INFO)
logging.disable(logging.WARNING)

import os.path
import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pychemkin.chemkin import *
from pychemkin.flux import *

job_pdd = ChemkinJob(
    name = 'PDD',
    chemFile = '',
    tempDir = '.',
)


# Temperatures 
#Tlist = [250, 325, 350, 375, 400, 425]  # C 
#
## Plotting styles
#colorwheel = ['#2c7fb8','#d7191c','#fdae61','#66c2a4']
##colorwheel = ['#253494','#2c7fb8','#41b6c4','#a1dab4']
##colorwheel = ['#4575b4','#66c2a4','#fec44f','#d7301f']
##colorwheel = ['#0868ac','#43a2ca','#7bccc4','#bae4bc']
#markers = ['o','s','^','D']
#lines = ['-','--','-.',':']
#
#    
#input = job_pdd.writeInputHomogeneousBatch(
#                      
#                      problemType = 'constrainPandT', # constrainVandSolveE, constrainPandSolveE, constrainPandT, constrainVandT                                   
#                      reactants=[(spc['PDD'],1),
#                                 ],     # Reactant (mole fraction)
#                      
#                      temperature = Tlist[0] + 273.15, # Temperature(K)
#                      pressure = 350,   # Pressure (bar)
#                      endTime = 5000*3600,   # End Time (sec), goes to 1000 hr
#                      
#                      Continuations = True,             # Continuations
#                      typeContinuation = 'NEWRUN',      # Type of continuation NEWRUN or CNTN
#                      Tlist=[T + 273.15 for T in Tlist[1:]] ,                  # Temperature (K) list of continuations
#                      rop = [spc['PDD']]
#                      #Plist=[350.0]#, 350.0, 350.0],                   # Pressure (atm) list of continuations
#                      
#                      #variableVolume = True,           ll # Variable volume true / false
#                      #variableVolumeProfile =vproFile,  # File name path for vpro (sec,cm3)
#                      
#                      #solverTimeStepProfile = solTimeStepFile # Solver time profile (sec)
#                      )
##if pdd:
##    job_pdd.preprocess()
##    job_pdd.run(input, model='CKReactorGenericClosed', pro=True)
##    job_pdd.postprocess(sens=False, rop=True, all=True, transpose=False)
#
## Extract Data
#
#time, molefracs = getMoleFraction(job_pdd.ckcsvFile, species=allSpecies)
#totalMoles = getTotalMoles(job_pdd.ckcsvFile)
#PDD_molefrac = molefracs[spc['PDD']]
from rmgpy.chemkin import loadChemkinFile

chemFile = 'chem.inp'
dictionary = 'species_dictionary.txt'
speciesList, reactionList = loadChemkinFile(chemFile, dictionary)


# Print off some ROP analysis
ROPData = {}
ROPData['PDD'] = extractROPData(job_pdd.ckcsvFile, species='PDD(1)', speciesList=speciesList, reactionList=reactionList)
TempIndex = 2
ROPData350 = ROPData['PDD'][TempIndex]
ROPData350.setFluxPercentages()
x = 10*3600   # 10 hours
ROPData350.sortTopRopReactions(x=x, numReactions=10)

ROPData350.sortTopRopSpecies(x=x, numSpecies=10)