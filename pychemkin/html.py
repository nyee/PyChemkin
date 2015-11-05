from rmgpy.chemkin import getSpeciesIdentifier


STYLE_HEADER = """
<style type="text/css">
        body {
            font-family: sans-serif;
        }
        a {
            color: #993333;
            text-decoration: none;
        }
        a:visited {
            color: #993333;
        }
        a:hover {
            text-decoration: underline;
        }
        table.speciesList, table.reactionList {
            width: 100%;
            border-collapse: collapse;
        }
        table.speciesList th, table.reactionList th {
            text-align: left;
            padding: 10px;
        }
        table.speciesList td, table.reactionList td {
            padding: 5px;
        }
        
        tr.reaction td, tr.reaction th{
            border-top: 1px solid #808080;
        }
        td.reactants {
            text-align: right;
        }
        td.products {
            text-align: left;
        }
        td.reactionArrow {
            text-align: center;
            font-size: 16px;
        }
        td.species img, td.reactants img, td.products img {
            vertical-align: middle;
        }
        tr.comment, td.comment {
           font-size: x-small;
           font-family: "Andale Mono", monospace;
        }
        tr.kinetics {
            font-size: small;
        }
        .KineticsData {
            # border: 1px solid gray;
        }
        .KineticsData th {
            width: 15em;
            word-wrap: none;
        }
        .KineticsData td {
            width: 3em;
        }

        .chemkin, .KineticsData_repr {
           white-space: pre-wrap;
           font-size: x-small;
           font-family: "Andale Mono", monospace;
        }

    </style>
    <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.4.1/jquery.min.js"></script>
    <script type="text/javascript" src="../../../external/jquery.min.js"></script>
   
"""

def renderSpecies(species, thermo=False):
    speciesLabel = getSpeciesIdentifier(species)
    speciesHTML = """
<tr class="species">
    <td class="label">{label}</td>
    <td class="structure"><a href={url}><img src="species/{filename}" alt="{label}" title="{label}"></a></td>
    <td class="smiles">{smiles}</td>
""".format(
           url=species.molecule[0].getURL(),
           filename=speciesLabel+'.png',
           label=speciesLabel,
           smiles=species.molecule[0].toSMILES(),
           )


    if thermo:
        thermo = species.thermo
        speciesHTML += """<td class="thermo">
        <table align="center">
            <tr>
                <th>H298</th>
                <th>S298</th>
                <th>Cp300</th>
                <th>Cp500</th>
                <th>Cp1000</th>
                <th>Cp1500</th>
            </tr>
            <tr>
                <td>{H298:.2f}</td>
                <td>{S298:.2f}</td>
                <td>{Cp300:.2f}</td>
                <td>{Cp500:.2f}</td>
                <td>{Cp1000:.2f}</td>
                <td>{Cp1500:.2f}</td>
            </tr>
            <tr class="comment"><th colspan=6>Comments</th></tr>
            <tr class="comment"><td colspan=6>{comments}</td></tr>
        </table></td>
""".format(H298=thermo.getEnthalpy(298)/4184,
           S298=thermo.getEntropy(298)/4.184,
           Cp300=thermo.getHeatCapacity(300)/4.184,
           Cp500=thermo.getHeatCapacity(300)/4.184,
           Cp1000=thermo.getHeatCapacity(300)/4.184,
           Cp1500=thermo.getHeatCapacity(300)/4.184, 
           comments=thermo.comment)

    speciesHTML += """</tr>"""
    
    return speciesHTML



def renderReactionROP(reactionROP, speciesList, showROP = True, idx=0):
    
    reaction = reactionROP.rxnObject
    
    reactionHTML = """
    <tr class="reaction"><th><a href="{url}" title="Search on RMG website" class="searchlink">{index}.</a> {rxnString}</th>
    """.format(rxnString=reaction.toChemkin(speciesList, kinetics=False), url=reaction.getURL(), index=reaction.index)
    
    if showROP:
        reactionHTML += """<td>Flux = {flux}</td>
        <td>Flux % = {percent}</td>""".format(flux=reactionROP.data[idx], percent=reactionROP.fluxPercentage[idx])
    reactionHTML += """</tr>"""
    
    reactantsHTML = []
    for reactant in reaction.reactants:
        reactantsHTML.append('<a href="{mol_url}"><img src="species/{label}.png" alt="{label}" title="{label}"></a>'.format(
                                                                                                                            mol_url=reactant.molecule[0].getURL(),
                                                                                                                            label=getSpeciesIdentifier(reactant),
                                                                                                                            ))
        
    reactionHTML += """<tr><td class="reactants">"""
    reactionHTML += ' + '.join(reactantsHTML)
    reactionHTML += """</td>"""
    
    reactionHTML += """<td class="reactionArrow">"""    
    if reaction.reversible:
        reactionHTML +="&hArr;"
    else:
        reactionHTML += "&rarr;"
    reactionHTML += """</td>"""
    
    productsHTML = []
    for product in reaction.products:
        productsHTML.append('<a href="{mol_url}"><img src="species/{label}.png" alt="{label}" title="{label}"></a>'.format(
                                                                                                                            mol_url=product.molecule[0].getURL(),
                                                                                                                            label=getSpeciesIdentifier(product),
                                                                                                                            ))
    
    reactionHTML += """<td class="products">"""
    reactionHTML += ' + '.join(productsHTML)
    reactionHTML += """</td></tr>"""
#    reactionHTML += """
#    <td class="family">{source}</td></tr>
#    """.format(source=reaction.getSource().label)
    
#    reactionHTML +="""
#    <tr class="kinetics">
#    <td colspan="4">{kinetics}</td>
#    </tr>""".format(kinetics=reaction.kinetics.toHTML())
    
    reactionHTML += """
    <tr class="chemkin">
    <th>CHEMKIN String</th></tr>
    <tr class="chemkin">
    <td>{chemkin}</td>
    </tr>
""".format(
           chemkin=reaction.toChemkin(speciesList)
           )

    return reactionHTML