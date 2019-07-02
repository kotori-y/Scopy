# -*- coding: utf-8 -*-
__doc__ = """
ScoPy(provisional), which based on rdkit and openbabel, 
has been designed for checking molecule(s) wheather or not 
has(have) some toxic substructure(s) through comparing the 
SMILES of molecule(s) and the SMARTS of toxic substructure(s) 
which obtained from ToxAlerts Database(https://ochem.eu/alerts/home.do).
========================================================================


Main Features
-------------
Here are just some slight works that we have done:
	
	-Checking molecule with different formats under the supporting of rdkit, 
	 including: SMILES, MolFile, SDMolSupplier.
	-If the SMILES cannot be recognized by rdkit, the SMILES will be converted to
	 Canonical SMILES through a function, based on openbabel module, 
	 which contributed by CatKin(https://blog.catkin.moe/).
	-Here you can check structural alerts for various toxicological endpoints
	 Endpoints including: 
	 0:Acute_Aquatic_Toxicity\t\t\t1:AlphaScreen-FHs
	 2:AlphaScreen-GST-FHs\t\t\t\t3:AlphaScreen-HIS-FHs
	 4:Biodegradable_compounds\t\t\t5:Chelating_agents
	 6:Custom_filters_\t\t\t\t7:Developmental_and_mitochondrial_toxicity
	 8:Extended_Functional_Groups_(EFG)\t\t9:Genotoxic_carcinogenicity,_mutagenicity
	 10:Idiosyncratic_toxicity_(RM_formation)\t11:LD50_mo_oral
	 12:Luciferase_Inhibitory_Activity\t\t13:Nonbiodegradable_compounds
	 14:Non-genotoxic_carcinogenicity\t\t15:PAINS_compounds
	 16:Potential_electrophilic_agents\t\t17:Promiscuity
	 18:Reactive,_unstable,_toxic\t\t\t19:Skin_sensitization
	 20:UNIFAC\t\t\t\t\t21:All
	 you could directly pass the index of endpoint as the parameter, the default parameter is "All".
	 -The result will be output in the format: "The XXX matched the SMARTS: YYY, endpoint: ZZZ",
	  the some function could show which part match the SMARTS, used in Ipython console.

Future works
-------------
Here some problems we do not resolve yet, and some works need we to finish:

	-Some SMARTS cannot be recognized.
	-Correcting the SMARTS.
	-Improving codes.
	-Giving a probability value of output.
"""
