import csv
import numpy

def extractROPData(ckcsvFile, species = '', speciesList=[], reactionList=[]):
    """
    Extract the ROP information from a chemkin run for a single species.
    Returns a list of SpeciesROP objects. Will only return a list containing a single SpeciesROP object
    if there are no continuations.
    
    Will also identify corresponding species and reaction objects from a loaded Chemkin file if 
    speciesList and reactionList are given.
    """
    from rmgpy.chemkin import getSpeciesIdentifier
    # Pre-convert the speciesList and reactionList to strings
    speciesStringList = []
    if speciesList:
        for spec in speciesList:
            speciesStringList.append(getSpeciesIdentifier(spec))
    
    reactionStringList = []
    if reactionList:
        for rxn in reactionList:
            rxnString = rxn.toChemkin(kinetics=False)
            if rxn.reversible:
                rxnString = rxnString.replace('=','<=>')
            reactionStringList.append(rxnString)
    
    timeData = {}; distanceData = {}; speciesROPData = {}
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        
        for row in reader:
            label = row[0].strip()
            tokens = label.split('_')
            if tokens[0] == 'Time':
                if 'Soln' in tokens[-1]:
                    tsolnNumber = int(tokens[-1].split('#')[-1])
                else:
                    tsolnNumber = 0
                tdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
                timeData[tsolnNumber] = tdata
                continue
            
            if tokens[0] == 'Distance':      
                if 'Soln' in tokens[-1]:
                    dsolnNumber = int(tokens[-1].split('#')[-1])
                else:
                    dsolnNumber = 0
                ddata = numpy.array([float(r) for r in row[2:]], numpy.float)
                dunits = row[1].strip()[1:-1].lower()
                ddata *= {'cm': 1.0, 'mm': 0.1, 'm': 100.}[dunits]
                distanceData[dsolnNumber] = ddata
                continue

            if len(tokens) > 1:
                if tokens[1] == 'ROP':
                    species_string = tokens[0]
                    if species == species_string:
                        # Create a SpeciesROP object for this solution number      
                        if 'Soln' in tokens[-1]:
                            solutionNumber = int(tokens[-1].split('#')[-1])
                        else:
                            # Only 1 solution set
                            solutionNumber = 0
                        if solutionNumber not in speciesROPData.keys():
                            speciesROPData[solutionNumber] = SpeciesROP(species, ropData=[])
                            if speciesList:
                                speciesROPData[solutionNumber].identifySpecies(speciesList, speciesStringList)
                        if tokens[3] == 'Total': 
                            # Total ROP rate
                            speciesROPData[solutionNumber].totalRopData = numpy.array([float(r) for r in row[2:]], numpy.float)
                        else:
                            rxnString = tokens[2]
                            rxnNumber = int(tokens[3].split('#')[-1])
                            rxnRopData = numpy.array([float(r) for r in row[2:]], numpy.float)
                            dataObject = ReactionROP(rxnString = rxnString, rxnNumber = rxnNumber, data = rxnRopData)
                            if reactionList:
                                dataObject.identifyReaction(reactionList, reactionStringList)
                            speciesROPData[solutionNumber].ropData.append(dataObject)
    
    speciesROPList = []
    if timeData:
        if len(timeData) != len(speciesROPData):
            raise Exception('Number of time and ROP data sets are mismatched.  Something went wrong during parsing')
        else:
            for solutionNumber in sorted(speciesROPData.keys()):
                speciesROP = speciesROPData[solutionNumber]
                speciesROP.xvarData = timeData[solutionNumber]
                speciesROPList.append(speciesROP)

    elif distanceData:
        if len(distanceData) != len(speciesROPData):
            raise Exception('Number of time and ROP data sets are mismatched.  Something went wrong during parsing')
        else:
            for solutionNumber in sorted(speciesROPData.keys()):
                speciesROP = speciesROPData[solutionNumber]
                speciesROP.xvarData = distanceData[solutionNumber]
                speciesROPList.append(speciesROP)

    else:
        raise Exception('Did not parse any time or distance data.  Something went wrong during parsing or simulation')
    
    return speciesROPList
    

class SpeciesROP(object):
    """
    A class for storing ROP data for a species at a single condition
    Stores the species information, the time data, total rop data, and rop data
    
    species should be a string
    species object should be an rmgpy Species object
    xvar data should be stored in the form of a numpy vector (could be either time or distance)
    total rop data should also be stored in the form of a numpy vector
    ropData comes as a list of ReactionROP objects, which each store rxn string, rxn number, and data
    """
    def __init__(self, species, xvarData=None, totalRopData=None, ropData=None):
        self.species = species
        self.xvarData = xvarData
        self.ropData = ropData
        self.totalRopData = totalRopData
        
        self.speciesObject = None  # To be set later
        self.totalPositiveRate = None
        self.totalNegativeRate = None
        
    def identifySpecies(self, speciesList, speciesStringList):
        """
        Returns a species object from a list of species based this object's species string
        Also stores it into the speciesObject
        """
        for i, spec in enumerate(speciesStringList):
            if self.species == spec:
                self.speciesObject = speciesList[i]
                return
        raise Exception('Did not find species {0} in species list.'.format(self.species))
            
    def getTotalRates(self):
        """
        Get the total positive and negative rates
        """
        self.totalPositiveRate = numpy.zeros_like(self.xvarData)
        self.totalNegativeRate = numpy.zeros_like(self.xvarData)
        for i in range(len(self.xvarData)):
            totalPositiveRate = 0
            totalNegativeRate = 0
            for rxnRop in self.ropData:
                rop = rxnRop.data[i]
                if rop <= 0:
                    totalNegativeRate += rop
                else:
                    totalPositiveRate += rop
            self.totalPositiveRate[i] = totalPositiveRate
            self.totalNegativeRate[i] = totalNegativeRate
        
    def setFluxPercentages(self):
        """
        Set the flux percentage for each reaction ROP data as a function of total positive or total negative rate
        """
        if not self.totalPositiveRate:
            self.getTotalRates()
        totalPositiveRate = self.totalPositiveRate
        totalNegativeRate = self.totalNegativeRate
        
        
        for rxnRop in self.ropData:
            fluxPercentage = numpy.zeros_like(self.xvarData)
            for i in range(len(rxnRop.data)):
                if rxnRop.data[i] <=0:
                    fluxPercentage[i] = rxnRop.data[i]/totalNegativeRate[i]*100
                else:
                    fluxPercentage[i] = rxnRop.data[i]/totalPositiveRate[i]*100
            rxnRop.fluxPercentage = fluxPercentage
            
    def sortTopRopReactions(self, x = 0, numReactions = 5, speciesList=[]):
        """
        Sort by the largest absolute flux percentages at a given value for the x variable (either time or distance).  
        The default number of reaction ROPs returned in either direction is set by numReactions
        """
        idx = (numpy.abs(numpy.array(self.xvarData) - x)).argmin()
        
        # Sort using the top absolute value of rop at that time or distance
        sortedROP = sorted(self.ropData, key=lambda rxnRop : abs(rxnRop.data[idx]), reverse=True)  

        sortedPositiveROP = []
        sortedNegativeROP = []
        for rxnRop in sortedROP:
            if rxnRop.data[idx] <= 0:
                sortedNegativeROP.append(rxnRop)
            else:
                sortedPositiveROP.append(rxnRop)
            
        print 'Actual value of x = {0}'.format(self.xvarData[idx])
        print 'Net ROP = {0}'.format(self.totalRopData[idx])        
        
        print '\nTotal Negative ROP = {0}'.format(self.totalNegativeRate[idx])                
        print '\nPrinting the top depletion reactions (negative flux) ...'
        
        for i, rxnRop in enumerate(sortedNegativeROP[:numReactions]):
            print '{num}. {rxnString}\t Actual flux = {flux}\t % of Total Negative ROP = {fluxpercent}'.format(num=i+1, rxnString=rxnRop.rxnString, flux=rxnRop.data[idx], fluxpercent=rxnRop.fluxPercentage[idx])
            
            
        print '\nTotal Positive ROP = {0}'.format(self.totalPositiveRate[idx])
        print '\nPrinting the top accumulation reactions (positive flux) ...'
        
        for i, rxnRop in enumerate(sortedPositiveROP[:numReactions]):
            print '{num}. {rxnString}\t Actual flux = {flux}\t % of Total Positive ROP = {fluxpercent}'.format(num=i+1, rxnString=rxnRop.rxnString, flux=rxnRop.data[idx], fluxpercent=rxnRop.fluxPercentage[idx])
        
        reactionList = [rxnRop.rxnObject for rxnRop in sortedNegativeROP[:numReactions]]
        saveFluxAnalysisHTML('fluxAnalysis.html', idx, self, negativeROPList = sortedNegativeROP[:numReactions], positiveROPList = sortedPositiveROP[:numReactions], speciesList=speciesList)
        
    def sortTopRopSpecies(self, x = 0, numSpecies = 5):
        """
        Sort species by the largest absolute flux percentages at a given value for the x variable (either time or distance).  
        The default number of species is set by numSpecies
        
        Note that this ignores species that are reacting on the same side as the parent species and participating with it in bimolecular
        reactions.  It only looks at species that produce our parent species are are produced by the parent species.
        """
        
        def addFlux(species, flux, rxnROP, positiveSpeciesROP, negativeSpeciesROP):
            if flux <= 0:
                if species in negativeSpeciesROP:
                    negativeSpeciesROP[species][0] += flux
                    negativeSpeciesROP[species][1].append(rxnROP)
                else:
                    negativeSpeciesROP[species] = [flux, [rxnROP]]
            else:                        
                if species in positiveSpeciesROP:
                    positiveSpeciesROP[species][0] += flux
                    positiveSpeciesROP[species][1].append(rxnROP)
                else:
                    positiveSpeciesROP[species] = [flux, [rxnROP]]
                    
        idx = (numpy.abs(numpy.array(self.xvarData) - x)).argmin()
        
        # We want to ignore species that are on the same 
        # Sort using the top absolute value of rop at that time or distance
        positiveSpeciesROP = {}
        negativeSpeciesROP = {}
        for rxnROP in self.ropData:
            if self.speciesObject in rxnROP.rxnObject.reactants:
                for spec in rxnROP.rxnObject.products:
                    # Change sign on flux since we are looking at products on other side
                    flux = -rxnROP.data[idx]
                    addFlux(spec, flux, rxnROP, positiveSpeciesROP, negativeSpeciesROP)
                #for spec in rxnROP.rxnObject.reactants:
                #    flux = rxnROP.data[idx]
                #   addFlux(spec, flux, positiveSpeciesROP, negativeSpeciesROP)
                        
                        
            elif self.speciesObject in rxnROP.rxnObject.products:
                for spec in rxnROP.rxnObject.reactants:
                    # Change sign on flux since we are looking at reactants on other side
                    flux = -rxnROP.data[idx]
                    addFlux(spec, flux, rxnROP, positiveSpeciesROP, negativeSpeciesROP)
                #for spec in rxnROP.rxnObject.products:
                #   flux = rxnROP.data[idx]
                #   addFlux(spec, flux, positiveSpeciesROP, negativeSpeciesROP)
                    
            else:
                raise Exception('Species {0} was not found in neither the reactants or products of this ROP reaction {1}.'.format(self.species, rxnROP.rxnString))
            
        # Sort based off absolute flux
        sortedNegativeSpecies = sorted(negativeSpeciesROP.items(), key=lambda x: abs(x[1][0]), reverse=True)
        sortedPositiveSpecies = sorted(positiveSpeciesROP.items(), key=lambda x: abs(x[1][0]), reverse=True)

        print 'Actual value of x = {0}'.format(self.xvarData[idx])
        #print 'Net ROP = {0}'.format(self.totalRopData[idx])        
        
        print '\nTotal Negative ROP = {0}'.format(self.totalNegativeRate[idx])
        print '\nPrinting the top species being generated (negative flux from parent species {0}) ...'.format(self.species)
                
        for i, speciesTuple in enumerate(sortedPositiveSpecies[:numSpecies]):
            print '{num}. {species}\t Actual flux = {flux}\t % of Total Negative ROP = {fluxpercent}'.format(num=i+1, species=speciesTuple[0], flux=speciesTuple[1][0], fluxpercent=-speciesTuple[1][0]/self.totalNegativeRate[idx]*100)
            sortedRxns = sorted(speciesTuple[1][1], key=lambda x: abs(x.data[idx]), reverse=True)  # Sort the rxns the species came from 
            for rxnRop in sortedRxns[:5]:
                print '\t{rxnString}\t Actual flux = {flux}\t % Contribution = {percent}'.format(rxnString=rxnRop.rxnString, flux=rxnRop.data[idx], percent=abs(rxnRop.data[idx]/speciesTuple[1][0])*100)
        print '\nTotal Positive ROP = {0}'.format(self.totalPositiveRate[idx])
        print '\nPrinting the top species being depleted (positive flux to form parent species {0}) ...'.format(self.species)

        for i, speciesTuple in enumerate(sortedNegativeSpecies[:numSpecies]):
            print '{num}. {species}\t Actual flux = {flux}\t % of Total Negative ROP = {fluxpercent}'.format(num=i+1, species=speciesTuple[0], flux=speciesTuple[1][0], fluxpercent=-speciesTuple[1][0]/self.totalPositiveRate[idx]*100)
            sortedRxns = sorted(speciesTuple[1][1], key=lambda x: abs(x.data[idx]), reverse=True)  # Sort the rxns the species came from 
            for rxnRop in sortedRxns[:5]:
                print '\t{rxnString}\t Actual flux = {flux}\t % Contribution = {percent}'.format(rxnString=rxnRop.rxnString, flux=rxnRop.data[idx], percent=abs(rxnRop.data[idx]/speciesTuple[1][0])*100)

        
class ReactionROP(object):
    """
    A class for storing ROP data
    """
    def __init__(self, rxnString='', rxnNumber=0, data=None):
        self.rxnString = rxnString
        self.rxnNumber = rxnNumber
        self.data = data
        
        self.rxnObject = None
        self.fluxPercentage = None 
        
    def identifyReaction(self, reactionList, reactionStringList):
        """
        Returns a reaction object from a list of reactions based this object's rxn string
        Also stores it into the rxnObject
        """
        for i, rxn in enumerate(reactionStringList):
            if self.rxnString == rxn:
                self.rxnObject = reactionList[i] 
                return
        raise Exception('Did not find reaction {0} in reaction list.'.format(self.rxnString))
    
class FluxDiagram(object):
    """
    A class for generating a flux diagram
    """
    def __init__(self, initialSpecies='', tolerance=0.1):
        self.initialSpecies = ''
        self.tolerance = tolerance
    
################################################################################

def saveFluxAnalysisHTML(path, idx, speciesROP, negativeROPList = [], positiveROPList = [], speciesList=[]):
    """
    Save the flux analysis for a given species to HTML.
    """
    import os
    from rmgpy.molecule.draw import MoleculeDrawer
    from rmgpy.chemkin import getSpeciesIdentifier
    from .html import renderSpecies, renderReactionROP, STYLE_HEADER
    
    try:
        import jinja2
    except ImportError:
        print "jinja2 package not found; HTML output will not be saved."
        return

    path = os.path.abspath(path)
    dirname = os.path.dirname(path)
    
            
    
    if not os.path.isdir(os.path.join(dirname,'species')):
        os.makedirs(os.path.join(dirname,'species'))

    for spec in speciesList:
        # Draw molecules 
        speciesLabel = getSpeciesIdentifier(spec)
        fstr = os.path.join(dirname, 'species', '{0}.png'.format(speciesLabel))
        if not os.path.exists(fstr):
            try:
                MoleculeDrawer().draw(spec.molecule[0], 'png', fstr)
            except IndexError:
                raise Exception("{0} species could not be drawn because it did not contain a molecular structure. Please recheck your files.".format(speciesLabel))
            
    #title = 'Flux Analysis for Species {label} at x = {x}'.format(label=getSpeciesIdentifier(species), x=speciesROP.xvarData[idx])
    #self.totalNegativeRate[idx]
    #print '{num}. {rxnString}\t Actual flux = {flux}\t % of Total Positive ROP = {fluxpercent}'.format(num=i+1, rxnString=rxnRop.rxnString, flux=rxnRop.data[idx], fluxpercent=rxnRop.fluxPercentage[idx])
    
    environment = jinja2.Environment()    
    environment.filters['renderSpecies'] = renderSpecies
    environment.filters['renderReactionROP'] = renderReactionROP
    environment.filters['getSpeciesIdentifier'] = getSpeciesIdentifier
    
    # Make HTML file
    template = environment.from_string(
"""<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN">
<html lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html;charset=utf-8" >
<title>Flux Analysis for Species {{speciesROP.speciesObject|getSpeciesIdentifier}} at x = {{speciesROP.xvarData[idx]}}</title>
{{style}}
</head>
<body>
<h1>Flux Analysis for Species {{speciesROP.speciesObject|getSpeciesIdentifier}} at x = {{speciesROP.xvarData[idx]}}</h1>
<h2>Species</h2>

<table class="speciesList">
<tr><th>Label</th><th>Structure</th><th>SMILES</th><th>Thermo</th>
{{speciesROP.speciesObject|renderSpecies(thermo=True)}}
</table>
<hr>

<h3>Net ROP = {{speciesROP.totalRopData[idx]}}</h3>

<hr size=3>
<h2>Top Negative Flux Reactions Depleting Species {{speciesROP.speciesObject|getSpeciesIdentifier}}</h2>
<h3>Total Negative ROP = {{speciesROP.totalNegativeRate[idx]}}</h3>

<table class="reactionList" hide_comment hide_kinetics hide_chemkin">
{% for rxnROP in negativeROPList %}
{{rxnROP|renderReactionROP(speciesList, showROP=True, idx=idx)}}
{% endfor %}
</table>

<hr size=3>
<h2>Top Positive Flux Reactions Accumulating Species {{speciesROP.speciesObject|getSpeciesIdentifier}}</h2>
<h3>Total Positive ROP = {{speciesROP.totalPositiveRate[idx]}}</h3>

<table class="reactionList" hide_comment hide_kinetics hide_chemkin">
{% for rxnROP in positiveROPList %}
{{rxnROP|renderReactionROP(speciesList, showROP=True, idx=idx)}}
{% endfor %}
</table>

</body>
</html>
""")

    f = open(path, 'w')
    f.write(template.render(style=STYLE_HEADER, 
                            idx = idx, 
                            speciesROP = speciesROP, 
                            negativeROPList = negativeROPList, 
                            positiveROPList = positiveROPList, 
                            speciesList=speciesList))
    f.close()
