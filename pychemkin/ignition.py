
################################################################################

def getIgnitionDelay(ckcsvFile, tol=1.0, species = []):
    """
    Return the ignition delay time from the given CKCSV file. This function
    uses (dP/dt)_max to find the ignition delay time. A ValueError is raised if
    this dP/dt value is below a certain threshold.
    
    Alternatively, provide a list of species which are maximized at tau.  
    If the list contains a single species label, ie. ['OH'] the ignition delay will be given as the
    time at which the concentraiton of that species is maximized.  If a list of species is given,
    then the multiplicative max of the concentrations of those species is returned.  ie. if 
    it is desired to calculated the ignition delay using max([CH][O]), then the species list 
    ['CH','O'] should be inputted.
    """
    
    tdata = None; Pdata = None; specdata = []
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            if label.startswith('Time'):
                tdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
            elif label.startswith('Pressure'):
                Pdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                Punits = row[1].strip()[1:-1].lower()
                Pdata *= {'dyne/cm2': 0.1, 'atm': 101325., 'Pa': 1.0, 'bar': 1e5, 'torr': 133.322368, 'mtorr': 0.133322368, 'psi': 6894.75729}[Punits] * 1e-5
            for spec in species:
                if spec in label:
                    specdata.append(numpy.array([float(r) for r in row[2:]], numpy.float))

    if tdata is None or Pdata is None:
        raise Exception('Unable to read time and/or pressure data from the given CKCSV file.')
    
    if species:
        if len(species) == 1:
            max_index = specdata[0].argmax()
            OHmetric = max(specdata[0])/2  # 0.0015136805  the particular value for GRI-Mech 3.0
            mindata = abs(specdata[0][0:max_index]-OHmetric)
            index = mindata.argmin()
            return tdata[index]
        else:
            
            multdata = numpy.ones(len(specdata[0]))
            ind1 = specdata[0].argmax()
            ind2 = specdata[1].argmax()
            for spec in specdata:  
                multdata *= spec
            print multdata
            print tdata
            print 'time max for species 1'
            print tdata[ind1]
            print 'time max for species 2'
            print tdata[ind2]
            index = multdata.argmax()
            print 'time max for multiplied together'
            print tdata[index]
            return tdata[index]
    else:
        dPdt = (Pdata[1:] - Pdata[:-1]) / (tdata[1:] - tdata[:-1])
        dPdt = dPdt[numpy.isfinite(dPdt)]
         
        #index = dPdt.argmax()
        index = next(i for i,d in enumerate(dPdt) if d==max(dPdt))
        if dPdt[index] < tol:
            raise ValueError('No ignition occurred in the given simulation.')
        
        return 0.5 * (tdata[index] + tdata[index+1])

def getIgnitionDelayOH(ckcsvFile, tol=1.0):
    """
    Return the ignition delay time from the given CKCSV file. This function
    uses (OH)_max to find the ignition delay time. 
    """
    
    tdata = None; OHdata = None
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            if label.startswith('Time'):
                tdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
            elif label.startswith('Mole_fraction_OH'):
                OHdata = numpy.array([float(r) for r in row[2:]], numpy.float)
    
    if tdata is None or OHdata is None:
        raise Exception('Unable to read time and/or OH data from the given CKCSV file.')
    
    OHdata = OHdata[numpy.isfinite(OHdata)] 
    index = OHdata.argmax()

    if index == len(OHdata)-1:
        raise ValueError('No ignition occurred in the given simulation.')

    return 0.5 * (tdata[index])

################################################################################

def getPeakOQOOHTime(ckcsvFile,spc, tol=1.0):
    """
    Return the ignition delay time from the given CKCSV file. This function
    uses (OQOOH)_max to find the ignition delay time. 
    """
    
    tdata = None; OQOOHdata = None
    spc_label =  'Mole_fraction_'+spc
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            if label.startswith('Time'):
                tdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
            elif label.startswith(spc_label):
                OQOOHdata = numpy.array([float(r) for r in row[2:]], numpy.float)
    
    if tdata is None or OQOOHdata is None:
        raise Exception('Unable to read time and/or OH data from the given CKCSV file.')
    
    OQOOHdata = OQOOHdata[numpy.isfinite(OQOOHdata)] 
    index = OQOOHdata.argmax()

    if index == len(OQOOHdata)-1:
        raise ValueError('No OQOOH peak found in the given simulation.')

    return (tdata[index])

################################################################################

def getIgnitionDelayStage1(ckcsvFile,stepsize = 1500,tol=1.0):
    """
    Return the ignition delay time from the given CKCSV file. This function
    uses (d2T/dt2) inflection to find the first stage ignition delay time. A ValueError is raised if
    this dT/dt value is below a certain threshold.
    """
    
    tdata = None; Tdata = None
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            if label.startswith('Time'):
                tdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
            elif label.startswith('Temperature'):
                Tdata = numpy.array([float(r) for r in row[2:]], numpy.float)
    
    if tdata is None or Tdata is None:
        raise Exception('Unable to read time and/or temperature data from the given CKCSV file.')
    
    from scipy import interpolate
    s = interpolate.interp1d(tdata[:], Tdata[:])
    xs = numpy.linspace(0.0, tdata[-1], stepsize)
    ys = s(xs)
    dT_dt = (ys[1:]-ys[:-1])/(xs[1:]-xs[:-1])
    xs1 = xs[1:]
    dT2_dt2 = (dT_dt[1:]-dT_dt[:-1])/(xs1[1:]-xs1[:-1])
    xs2 = xs1[1:]
    try:
        time_index_Tinflection = next(x[0] for x in enumerate(dT2_dt2) if x[1] < -100 )
    except StopIteration:
        raise ValueError('No T kink found.')
        #time_index_Tinflection = len(xs2)             
    
    Time_inflection = xs2[time_index_Tinflection]

    return Time_inflection
################################################################################


################################################################################

def runIgnitionReactionSensitivity(runChemkinJob, inputFile, dictionaryFile):
    """
    Supply a runChemkinJob python function which returns the ignition delay with a 
    chemkin file input.  This will run finite difference reaction sensitivities and 
    save them to a csv file.
    """
    
    from rmgpy.chemkin import loadChemkinFile, saveChemkinFile
    
    
    speciesList, reactionList = loadChemkinFile(inputFile, dictionaryPath = dictionaryFile, readComments = False)
    
    num_reactions = len(reactionList)
    
    factor_high = 1.05
    factor_low = 1. / factor_high
    
    worksheet = csv.writer(file('ignition_rxn_sensitivity.csv', 'w'))
    worksheet.writerow(['Index', 'Reaction', 'd[ln k]','tau_high','tau_low','d[ln tau]/d[ln k]'])
    
    print 'Running reaction sensitivity analysis using finite differences...'
    for index, reaction in enumerate(reactionList):
        rxn_index = index + 1
        rxn_string = reaction.toChemkin(kinetics = False)
        print 'At reaction {0} of {1}. {2}'.format(rxn_index, num_reactions, rxn_string)
                
        reaction.kinetics.changeRate(factor_high)
        saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose = False)
        tau_high = runChemkinJob('chem_temp.inp')
        reaction.kinetics.changeRate(1./factor_high)    # reset the kinetics
        
        reaction.kinetics.changeRate(factor_low)
        saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose = False)
        tau_low = runChemkinJob('chem_temp.inp')
        reaction.kinetics.changeRate(1./factor_low)     # reset the kinetics
        
        if tau_high != 0 and tau_low != 0:
            sens = numpy.log(tau_high / tau_low) / numpy.log(factor_high / factor_low)
        else:
            sens = 0
            
        worksheet.writerow([rxn_index, rxn_string, numpy.log(factor_high / factor_low), tau_high, tau_low, sens])
        
################################################################################

def runIgnitionThermoSensitivity(runChemkinJob, inputFile, dictionaryFile):
    """
    Supply a runChemkinJob python function which returns the ignition delay with a 
    chemkin file input.  This will run finite difference sensitivities to enthalpies and 
    save them to a csv file.
    """
    
    from rmgpy.chemkin import loadChemkinFile, saveChemkinFile, getSpeciesIdentifier
    from rmgpy.quantity import Quantity
    
    
    speciesList, reactionList = loadChemkinFile(inputFile, dictionaryPath = dictionaryFile, readComments = False)
    
    num_species = len(speciesList)
    
    deltaH = Quantity(0.5, 'kcal/mol').value_si
    
    worksheet = csv.writer(file('ignition_thermo_sensitivity.csv', 'w'))
    worksheet.writerow(['Species', 'd[del H] (kcal/mol)', 'tau_high', 'tau_low', 'd[ln tau]/d[del H]'])
    
    print 'Running thermo sensitivity analysis using finite differences...'
    for index, species in enumerate(speciesList):
        species_index = index + 1
        species_string = getSpeciesIdentifier(species)
        print 'At species {0} of {1}. {2}'.format(species_index, num_species, species_string)
        
        species.thermo.changeBaseEnthalpy(deltaH)
        saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose = False)
        tau_high = runChemkinJob('chem_temp.inp')
        species.thermo.changeBaseEnthalpy(-deltaH)    # reset the thermo
        
        species.thermo.changeBaseEnthalpy(-deltaH)
        saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose = False)
        tau_low = runChemkinJob('chem_temp.inp')
        species.thermo.changeBaseEnthalpy(deltaH)     # reset the kinetics
        
        if tau_high != 0 and tau_low != 0:
            sens = numpy.log(tau_high / tau_low) / (2 * deltaH)
        else:
            sens = 0
            
        worksheet.writerow([species_string, '1', tau_high, tau_low, sens])
