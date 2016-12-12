#!/usr/bin/env python
# encoding: utf-8

################################################
import logging
logging.disable(logging.INFO)
logging.disable(logging.WARNING)

import os.path
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import json
from pychemkin.chemkin import *
from rmgpy.chemkin import loadChemkinFile,getSpeciesIdentifier
from rmgpy.molecule import Molecule

#################################################

def getSpecListFromChemkinFiles(chemFile, dictionary):
  speciesList, reactionList = loadChemkinFile(chemFile, dictionary)
  for species in speciesList:
      species.generateResonanceIsomers()
  return speciesList


def getModelSpecDict(speciesList, speciesSMILES):
  """
  This method is trying to get a dict with keys of user-friendly 
  species names and with values of species labels or a list of 
  species labels in kinetic mechanism.
  """
  # Find the dictionary of corresponding model species
  spc_dict = {} 
  for spec_name in speciesSMILES.keys():
      SMILES = speciesSMILES[spec_name]
      if isinstance(SMILES,list):
          for SMILES_string in SMILES:
              mol = Molecule().fromSMILES(SMILES_string)
              for spec in speciesList:
                  if spec.isIsomorphic(mol):
                      if spec_name in spc_dict:
                          spc_dict[spec_name].append(getSpeciesIdentifier(spec))
                      else:
                          spc_dict[spec_name] = [getSpeciesIdentifier(spec)]
                      break
              # else:
                  # print 'Could not find species {0} in model'.format(spec_name)

      else:

          mol = Molecule().fromSMILES(SMILES)
          for spec in speciesList:
              if spec.isIsomorphic(mol):
                  spc_dict[spec_name] = getSpeciesIdentifier(spec)
                  break
          # else:
              # print 'Could not find species {0} in model'.format(spec_name)
  return spc_dict

  # set global settings
def init_plotting():
    plt.rcParams['figure.figsize'] = (8, 6)
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Helvetica'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = 0.8*plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.loc'] = 'best'
    plt.rcParams['axes.linewidth'] = 1

    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

init_plotting()

##############################################
# RERUN
pdd = True

# model
model = 'model_example'

resultPath = 'Results'
if not os.path.isdir(resultPath):
    os.makedirs(resultPath)

speciesSMILES = {'PDD':'CCCCCCCCCCCCc1ccccc1',
'toluene':'Cc1ccccc1',
'ethylbenzene':'CCc1ccccc1', 
'propylbenzene':'CCCc1ccccc1',
'butylbenzene':'CCCCc1ccccc1', 
'pentylbenzene':'CCCCCc1ccccc1',
'undecane':'CCCCCCCCCCC',
'undecene':['C=CCCCCCCCCC','CC=CCCCCCCCC','CCC=CCCCCCCC','CCCC=CCCCCCC','CCCCC=CCCCCC'],
'decane':['CCCCCCCCCC','C=CCCCCCCCC'],
'ethane':'CC',
'methane':'C',
'hexene':['CCCCC=C','CCCC=CC','CCC=CCC'],
'hexane':'CCCCCC',
'benzene':'c1ccccc1',
'heptene':['CCCCCC=C','CCCCC=CC','CCCC=CCC'],
'heptane':'CCCCCCC',
'octene':['CCCCCCC=C','CCCCCC=CC','CCCCC=CCC','CCCC=CCCC'],
'octane':'CCCCCCCC',
'nonane':'CCCCCCCCC',
'styrene':'C=CC1C=CC=CC=1',
}

chemFile = 'chem.inp'
dictionary = 'species_dictionary.txt'

print '\nSearching for species in {} model...'.format(model)
speciesList = getSpecListFromChemkinFiles(chemFile, dictionary)
spc = getModelSpecDict(speciesList, speciesSMILES)

allSpecies = []
for spec in speciesList:
    allSpecies.append(getSpeciesIdentifier(spec))


#######################
# Surrogate PYROLYSIS #
#######################

exptPath = os.path.abspath('./')
for exptDataFile in os.listdir(exptPath):
  if exptDataFile.endswith('.json'):
    try:
      split_tup = exptDataFile.split('.')[0].split('_')
      sample_id = int(split_tup[1][6:])
      Temp = int(split_tup[2][:-1])
    except:
      continue
  else:
    continue

  # condition path
  conditionPath = os.path.join(resultPath, 'new_sample{0}_{1}C'.format(sample_id, Temp))
  if not os.path.isdir(conditionPath):
      os.makedirs(conditionPath)

  expt_data_filename = exptDataFile

  print "Preparing simulation for new_sample{0} @ {1}C".format(sample_id, Temp)
  expt_mole_per_mass = json.load(open(expt_data_filename))

  # get initial mf
  initial_mole_per_mass = 0.0
  expt_init_mf = {}
  for spec_name in expt_mole_per_mass:
    if spec_name != 'Time':

      initial_mole_per_mass += expt_mole_per_mass[spec_name][0]

  for spec_name in expt_mole_per_mass:
    if spec_name != 'Time':

      expt_init_mf[spec_name] = expt_mole_per_mass[spec_name][0]/initial_mole_per_mass


  endTime = expt_mole_per_mass['Time'][-1] #hour

  job_pdd = ChemkinJob(
    name = 'new_surrogate',
    chemFile = chemFile,
    tempDir = model + '_MIT/' + "new_sample{0}_{1}C".format(sample_id, Temp),
  )

  input = job_pdd.writeInputHomogeneousBatch(
                        
                        problemType = 'constrainPandT', # constrainVandSolveE, constrainPandSolveE, constrainPandT, constrainVandT                                   
                        reactants=[(spc['PDD'],expt_init_mf['PDD']),
                        (spc['toluene'],expt_init_mf['toluene']),
                        (spc['undecane'],expt_init_mf['undecane']),
                                   ],     # Reactant (mole fraction)
                        
                        temperature = Temp + 273.15, # Temperature(K)
                        pressure = 300,   # Pressure (bar)
                        endTime = endTime*3600,   # End Time (sec), goes to 1000 hr
                        #Plist=[350.0]#, 350.0, 350.0],                   # Pressure (atm) list of continuations
                        
                        #variableVolume = True,           ll # Variable volume true / false
                        #variableVolumeProfile =vproFile,  # File name path for vpro (sec,cm3)
                        
                        #solverTimeStepProfile = solTimeStepFile # Solver time profile (sec)
                        )
  if pdd:
      job_pdd.preprocess()
      job_pdd.run(input, model='CKReactorGenericClosed', pro=True)
      job_pdd.postprocess(sens=False, rop=True, all=True, transpose=False)

  # Extract Data

  time, molefracs = getLumpedMoleFraction(job_pdd.ckcsvFile, spc)
  _, all_molefracs = getMoleFraction(job_pdd.ckcsvFile, allSpecies)
  totalMoles = getTotalMoles(job_pdd.ckcsvFile)
  PDD_molefrac = molefracs['PDD']

  # get simulation total mass
  simulation_total_mass = 0.0
  for spec_name in expt_init_mf:
    initial_totalMole = totalMoles[0][0]
    
    try:
      molWeight = Molecule().fromSMILES(speciesSMILES[spec_name]).getMolecularWeight()*1000
    except:
      molWeight = Molecule().fromSMILES(speciesSMILES[spec_name][0]).getMolecularWeight()*1000
    simulation_total_mass += molefracs[spec_name][0][0]*initial_totalMole*molWeight

  ###############################
  # Plot MIT Speciation Data #
  ###############################



  colorwheel = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e']
  markers = ['o','s','^','D','v']

  predicted_mfProfiles = {}
  tdata = numpy.array(time[0])
  idx = (numpy.abs(numpy.array(tdata)-endTime*3600.)).argmin()  # get time at 4 hrs

  
  #################
  ### Reactants ###
  #################
  plt.figure()
  ax = plt.gca()
  i=0
  legend=[]#'wo','wo']
  legendlabels=[]#'Experiment', 'Simulation']
  reactants_plot = ['PDD', 'toluene', 'undecane']
  for key in reactants_plot:
      trace, = plt.plot(expt_mole_per_mass['Time'], 1e3*numpy.array(expt_mole_per_mass[key]), markers[i], color=colorwheel[i], label=key + ' Experiment', markersize=10)  
      legend.append(trace)
      legendlabels.append(key)

      predicted_mfProfiles[key] = molefracs[key][0]*totalMoles[0]/simulation_total_mass
      trace, = plt.plot(tdata[:idx]/3600, 1e3*predicted_mfProfiles[key][:idx], label = key + ' Simulation', color=colorwheel[i], linewidth=3)
      legend.append(trace)
      legendlabels.append('')
      print key
      print "Experiment*1000: {}".format(expt_mole_per_mass[key][-1]*1000)
      print "Model*1000: {}".format(predicted_mfProfiles[key][-1]*1000)

      
      i = i+1
  plt.xlabel('Time (hr)')
  plt.ylabel('Product Moles Per Total Mass (moles/kg)')
  plt.xlim((0, endTime + 0.1*endTime))
  plt.title('sample{0}: Reactants at {1}$^\circ$C'.format(sample_id, Temp))
  ax.xaxis.labelpad = 10
  plt.tight_layout(h_pad=1.0, w_pad=2.0,rect=[0, 0, 1, 0.94])
  plt.figlegend(legend,legendlabels,loc='upper center', numpoints=1, bbox_to_anchor=(0.539, 0.87),  ncol=len(reactants_plot))
  plt.savefig(conditionPath + '/Reactants.png')


  ###############################
  ### Major Aromatic Products ###
  ###############################
  plt.figure()
  ax = plt.gca()
  i=0
  legend=[]#'wo','wo']
  legendlabels=[]#'Experiment', 'Simulation']
  aro_products_plot = ['styrene', 'ethylbenzene', 'propylbenzene']
  for key in aro_products_plot:
      trace, = plt.plot(expt_mole_per_mass['Time'], 1e3*numpy.array(expt_mole_per_mass[key]), markers[i], color=colorwheel[i], label=key + ' Experiment', markersize=10)  
      legend.append(trace)
      legendlabels.append(key)

      predicted_mfProfiles[key] = molefracs[key][0]*totalMoles[0]/simulation_total_mass
      trace, = plt.plot(tdata[:idx]/3600, 1e3*predicted_mfProfiles[key][:idx], label = key + ' Simulation', color=colorwheel[i], linewidth=3)
      legend.append(trace)
      legendlabels.append('')
      print key
      print "Experiment*1000: {}".format(expt_mole_per_mass[key][-1]*1000)
      print "Model*1000: {}".format(predicted_mfProfiles[key][-1]*1000)
      
      i = i+1
  plt.xlabel('Time (hr)')
  plt.ylabel('Product Moles Per Total Mass (moles/kg)')
  plt.xlim((0, endTime + 0.1*endTime))
  plt.title('sample{0}: Major Aromatic Products at {1}$^\circ$C'.format(sample_id, Temp))
  ax.xaxis.labelpad = 10
  plt.tight_layout(h_pad=1.0, w_pad=2.0,rect=[0, 0, 1, 0.94])
  plt.figlegend(legend,legendlabels,loc='upper center', numpoints=1, bbox_to_anchor=(0.539, 0.87),  ncol=len(aro_products_plot))
  plt.savefig(conditionPath + '/Major_AromaticProducts.png')
  
  ################################
  ### Major Aliphatic Products ###
  ################################
  plt.figure()
  ax = plt.gca()
  i=0
  legend=[]#'wo','wo']
  legendlabels=[]#'Experiment', 'Simulation']
  ali_products_plot = ['heptane', 'octane', 'nonane', 'decane']
  for key in ali_products_plot:
      trace, = plt.plot(expt_mole_per_mass['Time'], 1e3*numpy.array(expt_mole_per_mass[key]), markers[i], color=colorwheel[i], label=key + ' Experiment', markersize=10)  
      legend.append(trace)
      legendlabels.append(key)

      predicted_mfProfiles[key] = molefracs[key][0]*totalMoles[0]/simulation_total_mass
      trace, = plt.plot(tdata[:idx]/3600, 1e3*predicted_mfProfiles[key][:idx], label = key + ' Simulation', color=colorwheel[i], linewidth=3)
      legend.append(trace)
      legendlabels.append('')
      print key
      print "Experiment*1000: {}".format(expt_mole_per_mass[key][-1]*1000)
      print "Model*1000: {}".format(predicted_mfProfiles[key][-1]*1000)

      
      i = i+1
  plt.xlabel('Time (hr)')
  plt.ylabel('Product Moles Per Total Mass (moles/kg)')
  plt.xlim((0, endTime + 0.1*endTime))
  plt.title('sample{0}: Major Aliphatic Products at {1}$^\circ$C'.format(sample_id, Temp))
  ax.xaxis.labelpad = 10
  plt.tight_layout(h_pad=1.0, w_pad=2.0,rect=[0, 0, 1, 0.94])
  plt.figlegend(legend,legendlabels,loc='upper right', numpoints=1, bbox_to_anchor=(0.539, 0.87),  ncol=len(ali_products_plot)/2)
  plt.savefig(conditionPath + '/Major_AliphaticProducts.png')