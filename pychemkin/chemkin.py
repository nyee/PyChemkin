#!/usr/bin/env python
# encoding: utf-8

"""
This module contains several helper functions useful when running Chemkin jobs
from the command line.
"""

import os.path
import subprocess
import csv
import numpy
import sys

################################################################################

# The directory in which Chemkin is installed
chemkin_path = os.path.abspath(os.path.join(os.path.dirname(__file__),'chemkin_path'))

CHEMKIN_DIR = ''
if os.path.exists(chemkin_path):
    with open(chemkin_path, 'r') as f:
        path = f.readline()
        if not path:
            print "Warning: CHEMKIN_DIR has not been set properly.  Please do so in the text file 'chemkin_path'"
        else:
            CHEMKIN_DIR = os.path.abspath(path)
    
# The preamble to each Chemkin execution shell script
CHEMKIN_SCRIPT_PREAMBLE = """#!/bin/sh -v

# Define Chemkin running environment
. {0}

""".format(os.path.join(CHEMKIN_DIR, 'bin', 'chemkinpro_setup.ksh'))


class ChemkinJob(object):
    """
    A single (gas-phase) Chemkin job.
    """
    
    def __init__(self, name, chemFile, tempDir):
        self.name = name
        self.chemFile = os.path.abspath(chemFile)
        self.tempDir = os.path.abspath(tempDir)
        if not os.path.exists(self.tempDir):
            os.makedirs(self.tempDir)
        
    @property
    def ascFile(self):
        return os.path.join(self.tempDir, '{0}_gas.asc'.format(self.name))
    
    @property
    def inputFile(self):
        return os.path.join(self.tempDir, '{0}.inp'.format(self.name))
    
    @property
    def outputFile(self):
        return os.path.join(self.tempDir, '{0}.out'.format(self.name))
    
    @property
    def monFile(self):
        return os.path.join('{0}.out.mon'.format(self.name))
    
    @property
    def ckcsvFile(self):
        return os.path.join(self.tempDir, 'CKSoln.ckcsv')

    @property
    def preferenceFile(self):
        return os.path.join(self.tempDir, 'CKSolnList.txt')
    
    @property
    def dataZipFile(self):
        return os.path.join(self.tempDir, 'XMLdata_{0}.zip'.format(self.name))    
    
    def preprocess(self):
        """
        Run the Chemkin preprocessor on the chemistry input file.
        """
        # Write the preprocessor execution script to a file
        scriptFile = os.path.join(self.tempDir, 'RUNJOB_CKPreProcess.sh')
        with open(scriptFile, 'w') as stream:
            stream.write(CHEMKIN_SCRIPT_PREAMBLE)
            stream.write('# Run the Chemkin preprocessor\n')
            stream.write('{0} -i {1} -o {2}_gas.out -c {2}_gas.asc\n\n'.format(
                os.path.join(CHEMKIN_DIR, 'bin', 'chem'),
                os.path.abspath(self.chemFile),
                os.path.join(self.tempDir, self.name),
            ))
            
        # Execute the preprocessor script
        process = subprocess.Popen(('/bin/sh', scriptFile), cwd=self.tempDir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            print stdout
            print stderr
            quit()
        
    def run(self, jobParams, model, pro=True):
        """
        Run a Chemkin job.
        """
        dtdFileSrc = os.path.join(CHEMKIN_DIR, 'data', 'chemkindata.dtd')
        dtdFileDest = os.path.join(self.tempDir, 'chemkindata.dtd')
        scriptFile = os.path.join(self.tempDir, 'RUNJOB_CKRunProcessor.sh')
        # Write the job parameters to a file
        with open(self.inputFile, 'w') as stream:
            stream.write(jobParams)
        # Write the job execution script to a file
        with open(scriptFile, 'w') as stream:
            stream.write(CHEMKIN_SCRIPT_PREAMBLE)
            stream.write('# Delete any intermediate files left over from a previous run\n')
            stream.write('rm -f {0}\n\n'.format(self.dataZipFile))
            stream.write('rm -f {0}\n\n'.format(dtdFileDest))
            stream.write('# Copy files to working directory\n')
            stream.write('ln -s {0} {1}\n'.format(dtdFileSrc, dtdFileDest))
            stream.write('# Run the job\n')
            if pro:
                stream.write('CHEMKIN_MODE=Pro\n')
                stream.write('export CHEMKIN_MODE\n')
                stream.write('{0} -i {1} -o {2} -x {3} Pro -m {4} -c {5}\n'.format(model, self.inputFile, self.outputFile, self.dataZipFile, self.monFile, self.ascFile))
            else:
                stream.write('{0} -i {1} -o {2} -x {3} -m {4} -c {5}\n'.format(model, self.inputFile, self.outputFile, self.dataZipFile, self.monFile, self.ascFile))
        # Execute the job script
        subprocess.call(('/bin/sh', scriptFile), cwd=self.tempDir)

    def postprocess(self, sens=False, rop=False, all=False, transpose=True):
        """
        Run the Chemkin postprocessor.
        """
        # Set default preferences and write them to the preference file
        with open(self.preferenceFile, 'w') as stream:
            stream.write('VARIABLE VAR ALL\n')
            stream.write('VARIABLE SEN ALL\n')
            stream.write('VARIABLE ROP ALL\n') 
            stream.write('VARIABLE volume  1  0  0\n')
            stream.write('VARIABLE external_surface_area  0  0  0\n')
            stream.write('VARIABLE volumetric_heat_production_rate  0  0  0\n')
            stream.write('VARIABLE surface_temperature  0  0  0\n')
            stream.write('VARIABLE heat_loss_rate  0  0  0\n')
            stream.write('VARIABLE temperature  1  0  0\n')
            stream.write('VARIABLE mass  0  0  0\n')
            stream.write('VARIABLE pressure  1  0  0\n')
            stream.write('VARIABLE molar_conversion  0  0  0\n')
            stream.write('VARIABLE net_reaction_rate  0  0  0\n')
            stream.write('VARIABLE heat_production_per_reaction  0  0  0\n')
            stream.write('VARIABLE molecular_weight  0  0  0\n')
            stream.write('VARIABLE mass_density  0  0  0\n')
            stream.write('VARIABLE mixture_enthalpy  0  0  0\n')
            stream.write('VARIABLE sensible_enthalpy  0  0  0\n')
            stream.write('VARIABLE unburned_hydrocarbons  0  0  0\n')
            stream.write('VARIABLE volatile_organic_compounds  0  0  0\n')
            stream.write('VARIABLE exhaust_in_ppmvd  0  0  0\n')
            stream.write('VARIABLE all_single_point_values  0  0  0\n')
            stream.write('VARIABLE use_reaction_string_instead_of_index  1  0  0\n')
            stream.write('UNIT  Time  (sec)\n')
            stream.write('UNIT  Temperature  (K)\n')
            stream.write('UNIT  Pressure  (bar)\n')
            stream.write('UNIT  Volume  (cm3)\n')

        # Write the postprocessor execution script to a file
        scriptFile = os.path.join(self.tempDir, 'RUNJOB_CKPostProcessor.sh')
        with open(scriptFile, 'w') as stream:
            stream.write(CHEMKIN_SCRIPT_PREAMBLE)
            stream.write('# Extract solution data to CKCSV\n')
            stream.write('GetSolution {1}{2}{3}{4}{0}\n\n'.format(
                self.dataZipFile, 
                '-nosen ' if not sens else '',
                '-norop ' if not rop else '',
                '-all ' if all else '',
                '-p'+' CKSolnList.txt ',
            ))
            if transpose:
                stream.write('# Transpose the data to CSV\n')
                stream.write('CKSolnTranspose -column 500 {0}\n\n'.format(self.ckcsvFile))
            stream.write('# Delete postprocessor log file\n')
            stream.write('rm -f {0}\n\n'.format(os.path.join(self.tempDir, 'ckpp_*.log')))
        
        # Execute the postprocessor script
        subprocess.call(('/bin/sh', scriptFile), cwd=self.tempDir)
        
    def writeInputHomogeneousBatch(self,problemType, reactants, temperature, pressure, endTime, 
                      Continuations=False, typeContinuation = None, Tlist = [], Plist = [],
                      variableVolume=False, variableVolumeProfile = None, 
                      solverTimeStepProfile = None, sensitivity=[], rop=[]):
        """
        Write input file for homogeneous batch reactor
        """
        
        # Problem definition block 
                    
        input_stream=("""
! 
! problem type definition
!""")

        if problemType.lower() == 'constrainVandSolveE'.lower():
            input_stream+=("""
CONV   ! Constrain Volume And Solve Energy Equation
ENRG   ! Solve Gas Energy Equation""")
        elif problemType.lower() == 'constrainPandSolveE'.lower():
            input_stream+=("""
ENRG   ! Solve Gas Energy Equation""")
        elif problemType.lower() == 'constrainPandT'.lower():
            input_stream+=("""
TGIV   ! Constrain Pressure And Temperature
""")
        elif problemType.lower() == 'constrainVandT'.lower():
            input_stream+=("""
COTV   ! Constrain Volume And Temperature
TGIV   ! Fix Gas Temperature""")
             
        # Solver type definition block 
        
        input_stream+=("""
TRAN   ! Transient Solver""")
        
        # Physical property
        
        input_stream+=("""
! 
! physical property
! 
!Surface_Temperature   ! Surface Temperature Same as Gas Temperature
IFAC 0.1   ! Ignition Noise Filtering Factor
""")
            
        input_stream+=('PRES {0:g}   ! Pressure (atm)\n'.format(pressure/1.01325))
        input_stream+=('QLOS 0.0   ! Heat Loss (cal/sec)\n')
        input_stream+=('TEMP {0:g}   ! Temperature (K)'.format(temperature))
        
        if variableVolume:
            with open(variableVolumeProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    time = float(row[0].split(',')[0]) # (sec)
                    vol = float(row[0].split(',')[1]) # (cm3)
                    input_stream+=("""
VPRO {0:g} {1:g}   ! Volume (cm3)""".format(time,vol))
                       
        
        # Species property block
        
        input_stream+=("""
!
! species property
!
""")
        
        for reac , conc in reactants:            
            input_stream+=('REAC {0} {1:g} ! Reactant Fraction (mole fraction) \n'.format(reac,conc))
            
        # Solver control block
        
        input_stream+=("""! 
! solver control
! 
ADAP   ! Save Additional Adaptive Points
ASTEPS 20   ! Use Solver Integration Steps
ATLS 1.0E-6   ! Sensitivity Absolute Tolerance
ATOL 1.0E-20   ! Absolute Tolerance
MAXIT 4   ! Maximum Number of Iterations
RTLS 0.0001   ! Sensitivity Relative Tolerance
RTOL 1.0E-8   ! Relative Tolerance""")
        
        if solverTimeStepProfile:
            with open(solverTimeStepProfile, 'rb') as csvfile2:
                timeStepReader = csv.reader(csvfile2, delimiter=' ', quotechar='|')
                for row in timeStepReader:
                    time = float(row[0].split(',')[0]) # (sec)
                    vol = float(row[0].split(',')[1]) # (sec)
                    input_stream+=("""
STPTPRO {0:g} {1:g}           ! Solver Maximum Step Time (sec)""".format(time,vol))
            
        
        input_stream+=("""
TIME {0:g}                 ! End Time (sec)""".format(endTime))
        
        input_stream+=("""
! 
! output control and other misc. property
!
""")

        for species in sensitivity:
            input_stream+=('ASEN {0}    ! A-factor Sensitivity\n'.format(species)) 

        input_stream+=("""
EPSR 0.01   ! Threshold for Rate of Production
EPSS 0.001   ! Threshold for Species Sensitivity
EPST 0.001   ! Threshold for Temperature Sensitivity
GFAC 1.0   ! Gas Reaction Rate Multiplier
PRNT 1   ! Print Level Control
SIZE 10000000   ! Solution Data Block Size (bytes)
""")

        for species in rop:
            input_stream+=('ROP {0}    ! Rate of Production\n'.format(species))
              
        if Continuations:
            if numpy.array(Tlist).size:                
                for i in range(numpy.array(Tlist).shape[0]):
                    input_stream+=("""
{0}
END
TEMP {1:g}""".format(typeContinuation,numpy.array(Tlist)[i]))
            
            if numpy.array(Plist).size:         
                for i in range(numpy.array(Plist).shape[0]):
                    input_stream+=("""
{0}
END
PRES {1:g}""".format(typeContinuation,numpy.array(Plist)[i]/1.01325))
                                            
        input_stream+=('\nEND')
        
        return input_stream      

    def writeInputPlugFlow(self,problemType, reactants, 
                           startingAxialPosition, endingAxialPosition, diameter, momentum=True, massflowrate = None, sccmflowrate = None, 
                           temperature  = None, pressure  = None,  
                           temperatureProfile = None, pressureProfile = None, sensitivity=[], rop=[]):
    
        """
        Write input file for typical plug flow
        """
        
        # Problem definition block 
                    
        input_stream=("""
! 
! problem type definition
!
""")

        if momentum:
            input_stream+=("""
MOMEN ON   ! Turn on Momentum Equation
""")
        else:
            input_stream+=("""
MOMEN OFF   ! Turn off Momentum Equation
""")

        input_stream+=("""
PLUG   ! Plug Flow Reactor
RTIME ON   ! Turn on Residence Time Calculation
""")

        if problemType.lower() == 'FixGasTemperature'.lower():
            input_stream+=("""
TGIV   ! Fix Gas Temperature""")
        elif problemType.lower() == 'solveGasEnergy'.lower():
            input_stream+=("""
ENRG   ! Solve Gas Energy Equation""")
    
        
        # Physical property
        
        input_stream+=("""
! 
! physical property
! 
!Surface_Temperature   ! Surface Temperature Same as Gas Temperature
""")
        if massflowrate:
            input_stream+=('FLRT {0:g}   ! Mass Flow Rate (g/sec)\n'.format(massflowrate))
        elif sccmflowrate:
            input_stream+=('SCCM {0:g}   ! Volumetric Flow Rate in SCCM (standard-cm3/min@298.15K)\n'.format(sccmflowrate))
        else:
            raise Exception('Must indicate either a mass or sccm flow rate')
        
        if pressureProfile:
            with open(pressureProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    pos = float(row[0].split(',')[0]) # (cm)
                    pres = float(row[0].split(',')[1])*1.0/1.01325 # (atm)
                    input_stream+=("""
PPRO {0:g} {1:g}   ! Pressure (atm)""".format(pos,pres))
        else:
            input_stream+=('PRES {0:g}   ! Pressure (atm)\n'.format(pressure/1.01325))
            
        if temperatureProfile:
            with open(temperatureProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    pos = float(row[0].split(',')[0]) # (cm)
                    temp = float(row[0].split(',')[1]) # (K)
                    input_stream+=("""
TPRO {0:g} {1:g}   ! Temperature (K)""".format(pos,temp))
        else:
            input_stream+=('TEMP {0:g}   ! Temperature (K)'.format(temperature)) 
        
        # reactor dimension definition
        
        input_stream+=("""
VIS 0.0   ! Mixture Viscosity at Inlet (g/cm-sec)
! 
! reactor dimension definition
! 
""")
        input_stream+=('DIAM {0:g}   ! Diameter (cm)\n'.format(diameter))
        input_stream+=('XEND {0:g}   ! Ending Axial Position (cm)\n'.format(endingAxialPosition))
        input_stream+=('XSTR {0:g}   ! Starting Axial Position (cm)\n'.format(startingAxialPosition))
                       
        
        # Species property block
        
        input_stream+=("""
!
! species property
!
""")
        
        for reac , conc in reactants:            
            input_stream+=('REAC {0} {1:g} ! Reactant Fraction (mole fraction) \n'.format(reac,conc))
            
        # Solver control block
        
        input_stream+=("""! 
! 
! solver control
! 
ACHG 0.0   ! Maximum Absolute Change in Site Fractions
ADAP   ! Save Additional Adaptive Points
ATLS 1.0E-6   ! Sensitivity Absolute Tolerance
ATOL 1.0E-9   ! Absolute Tolerance
MAXIT 4   ! Maximum Number of Iterations
NNEG   ! Force Non-negative Solution
PSV 1.0E-8   ! Scaling Factor for Relaxing Surface Equations (cm/sec)
RCHG 1.0E-6   ! Maximum Relative Change in Site Fractions
RTLS 0.0001   ! Sensitivity Relative Tolerance
RTOL 1.0E-6   ! Relative Tolerance
TSTP 1.0   ! Initial Integration Step (cm""")          

        input_stream+=("""
! 
! output control and other misc. property
!
""")
        for species in sensitivity:
            input_stream+=('ASEN {0}    ! A-factor Sensitivity\n'.format(species)) 

        input_stream+=("""
EPSR 0.01   ! Threshold for Rate of Production
EPSS 0.001   ! Threshold for Species Sensitivity
EPST 0.001   ! Threshold for Temperature Sensitivity
GFAC 1.0   ! Gas Reaction Rate Multiplier
PRNT 1   ! Print Level Control
SIZE 10000000   ! Solution Data Block Size (bytes)
""")

        for species in rop:
            input_stream+=('ROP {0}    ! Rate of Production\n'.format(species))
                                            
        input_stream+=('\nEND')
        
        return input_stream      
    
    def writeInputJSR(self,problemType, reactants, tau,endtime, volume, 
                           temperature  = None, pressure  = None,
                           Continuations=False, typeContinuation = None, Tlist = [], Plist = [],                           
                           temperatureProfile = None, pressureProfile = None, sensitivity=[], rop=[]):
    
        """
        Write input file for JSR
        """
        
        # Problem definition block 
                    
        input_stream=("""
! 
! problem type definition
!
TRAN   ! Transient Solver 
""")

        if problemType.lower() == 'FixGasTemperature'.lower():
            input_stream+=("""
TGIV   ! Fix Gas Temperature""")
        elif problemType.lower() == 'solveGasEnergy'.lower():
            input_stream+=("""
ENRG   ! Solve Gas Energy Equation""")
    
        
        # Physical property
        
        input_stream+=("""
! 
! physical property
! 
!Surface_Temperature   ! Surface Temperature Same as Gas Temperature
""")
                
        if pressureProfile:
            with open(pressureProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    pos = float(row[0].split(',')[0]) # (cm)
                    pres = float(row[0].split(',')[1])*1.0/1.01325 # (atm)
                    input_stream+=("""
PPRO {0:g} {1:g}   ! Pressure (atm)""".format(pos, pres))
        else:
            input_stream+=('PRES {0:g}   ! Pressure (atm)\n'.format(pressure/1.01325))
        
        # Residence time 
        input_stream+=('TAU {0:g}   ! Residence time (sec)\n'.format(tau))

        if temperatureProfile:
            with open(temperatureProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    pos = float(row[0].split(',')[0]) # (cm)
                    temp = float(row[0].split(',')[1]) # (K)
                    input_stream+=("""
TPRO {0:g} {1:g}   ! Temperature (K)""".format(pos,temp))
        else:
            input_stream+=('TEMP {0:g}   ! Temperature (K)\n'.format(temperature)) 
        
        # volume
        input_stream+=('VOL {0:g}   ! Volume (cm3)\n'.format(volume))
        
        # Species property block
        
        input_stream+=("""
!
! species property
!
""")
        
        for reac , conc in reactants:            
            input_stream+=('REAC {0} {1:g} ! Reactant Fraction (mole fraction) \n'.format(reac,conc))
        # For transient solver you also need estimates of species initial gas fraction
        for reac , conc in reactants:            
            input_stream+=('XEST {0} {1:g} ! Initial Gas Fraction (mole fraction) \n'.format(reac,conc))
        
        # solver control    
        input_stream+=("""
! 
! solver control
! 
ADAP   ! Save Additional Adaptive Points
""")    
        input_stream+=('TIME {0:g}   ! End Time (sec) \n'.format(endtime))    
    
        input_stream+=("""
! 
! output control and other misc. property
!
""")

        for species in sensitivity:
            input_stream+=('ASEN {0}    ! A-factor Sensitivity\n'.format(species)) 

        input_stream+=("""
EPSR 0.01   ! Threshold for Rate of Production
EPSS 0.001   ! Threshold for Species Sensitivity
EPST 0.001   ! Threshold for Temperature Sensitivity
GFAC 1.0   ! Gas Reaction Rate Multiplier
PRNT 1   ! Print Level Control
SIZE 10000000   ! Solution Data Block Size (bytes)
""")

        for species in rop:
            input_stream+=('ROP {0}    ! Rate of Production\n'.format(species))
                                            
        if Continuations:
            if numpy.array(Tlist).size:                
                for i in range(numpy.array(Tlist).shape[0]):
                    input_stream+=("""
{0}
END
TEMP {1:g}""".format(typeContinuation,numpy.array(Tlist)[i]))
            
            if numpy.array(Plist).size:         
                for i in range(numpy.array(Plist).shape[0]):
                    input_stream+=("""
{0}
END
PRES {1:g}""".format(typeContinuation,numpy.array(Plist)[i]/1.01325))
                                            
        input_stream+=('\nEND')        
        return input_stream      


################################################################################

def getVariable(ckcsvFile, variable=None):
    """
    Return the time and variable data (i.e. Pressure, Volume, or Temperature) from the given CKCSV file.  Returns the data
    in the form of [time_soln1, timesoln2, ...] [variable_data_soln1, variable_data_soln2,...]
    Returns time in units of [seconds]. Returns temperature in units of [K].
    Returns pressure in units of [bar]. Returns volume in units of [cm3].
    """
    if variable.lower() not in ['pressure','volume','temperature']:
        raise Exception('Can only parse Pressure, Volume, or Temperature from CKCSV file.')

    timeData = []; varData = []

    units = {'k': 1.0, 'bar': 1, 'mpa': 10, 'pa': 1e-5, 'atm': 1.01325, 'cm3': 1.0, 'm3': 1e6}
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            tokens = label.split('_')
            if tokens[0] == 'Time':
                tdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
                timeData.append(tdata)      

            if tokens[0].lower() == variable.lower():
                vdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                vunits = row[1].strip()[1:-1].lower()
                vdata = vdata*units[vunits]
                varData.append(vdata)
    
    return timeData, varData

def getMoleFraction(ckcsvFile, species=[]):
    """
    Return the time and mole fraction data from the given CKCSV file.  Returns the data
    in the form of [time_soln1, timesoln2, ...] specdata[species_name]=[spec_molefrac_soln1, spec_molefrac_soln2,...].
    Time is returned in units of [seconds].
    """

    timeData = []; distanceData = []; specData = {}

    for spec in species:
        specData[spec] = []
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            tokens = label.split('_')
            if tokens[0] == 'Time':
                tdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
                timeData.append(tdata)
            
            if tokens[0] == 'Distance':      
                ddata = numpy.array([float(r) for r in row[2:]], numpy.float)
                dunits = row[1].strip()[1:-1].lower()
                ddata *= {'cm': 1.0, 'mm': 0.1, 'm': 100.}[dunits]
                distanceData.append(ddata)

            if label.startswith('Mole_fraction'):
                for spec in species:
                    if tokens[2] == spec:
                        specData[spec].append(numpy.array([float(r) for r in row[2:]], numpy.float))
    if timeData:
        return timeData, specData
    elif distanceData:
        return distanceData, specData
    else:
        return

def getTotalMoles(ckcsvFile):
    """
    Return the time and total moles data from the given CKCSV file.  
    Returns the data in the form of [time_soln1, timesoln2, ...], [total_moles_soln1, total_moles_soln2, ...]
    Moles is returned in units of [moles].
    """
    
    tempData = []; presData = []; volData = []; totalMolesData = []

    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            tokens = label.split('_')

            if  tokens[0] == 'Temperature':
                Tdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                Tunits = row[1].strip()[1:-1].lower()
                try:
                    Tdata *= {'k': 1.0}[Tunits] 
                except KeyError:
                    print 'Unknown units {0} for temperature. Cannot obtain total moles.'.format(Tunits)
                tempData.append(Tdata)

            if tokens[0] == 'Pressure': 
                Pdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                Punits = row[1].strip()[1:-1].lower()
                try:
                    Pdata *= {'bar': 0.1, 'mpa': 1.0, 'pa': 1e-6, 'atm': 0.101325}[Punits]
                except KeyError:
                    print 'Unknown units {0} for pressure. Cannot obtain total moles.'.format(Punits)
                presData.append(Pdata)
            
            if tokens[0] == 'Volume':                
                Vdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                Vunits = row[1].strip()[1:-1].lower()
                try:
                    Vdata *= {'cm3': 1.0, 'm3': 1e6}[Vunits]
                except KeyError:
                    print 'Unknown units {0} for volume. Cannot obtain total moles.'.format(Vunits)
                volData.append(Vdata)

    R = 8.3145
    for i in range(len(tempData)):
        totalMolesData.append(presData[i]*volData[i]/R/tempData[i])

    return totalMolesData
    