#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Internet connection required

import pubchempy as pcp
import ipywidgets as widgets
from IPython.display import Image
import requests

# Textbox for Compound search
textbox_compound = widgets.Text(
    placeholder='Enter Compound (e.g. Aspirin)',
    disabled=False   
)
display(textbox_compound)


# In[2]:


# First button
button_getcompound = widgets.Button(
    description='Get Compound',
    disabled=False
)
# Output widgets for displaying the output
output = widgets.Output()
output2 = widgets.Output()
display(button_getcompound)

# Create Textboxes here, assign values and display in button function
textbox_iupac = widgets.Text(
    description='IUPAC name:',
    disabled=False   
)

textbox_formula = widgets.Text(
    description='Molecular formula:',
    disabled=False   
)

textbox_weight = widgets.Text(
    description='Molecular weight:',
    disabled=False   
)

textbox_xlogp = widgets.Text(
    description='XlogP:',
    disabled=False   
)

textbox_smiles = widgets.Text(
    description='SMILES:',
    disabled=False   
)

textbox_db = widgets.Text(
    description='Double bonds:',
    disabled=False   
)

textbox_tb = widgets.Text(
    description='Triple bonds:',
    disabled=False   
)

textbox_rings = widgets.Text(
    description='Rings:',
    disabled=False   
)

# Search the compound from textbox in pubchem database and display chemical properties
def get_compound_from_string(b):
    output.clear_output()
    output2.clear_output()
    with output:
        try:
            compound_list = pcp.get_compounds(textbox_compound.value, 'name')
            c = compound_list[0]
        except:
            eingabe = "\"" + textbox_compound.value + "\""
            print("Compound", eingabe, "not found. Please try again.")
        else:
            # Count double and triple bonds
            doublebonds = c.isomeric_smiles.count('=') 
            triplebonds = c.isomeric_smiles.count('#')
            
            # Count rings - multiple connected rings only count as one ring
            rings = int(c.isomeric_smiles.count('1') / 2)

            # Set values and display textboxes
            textbox_iupac.value = c.iupac_name
            display(textbox_iupac)
            
            textbox_formula.value = c.molecular_formula
            display(textbox_formula)
            
            textbox_weight.value = str(c.molecular_weight)
            display(textbox_weight)
            
            textbox_xlogp.value = str(c.xlogp)
            display(textbox_xlogp)
        
            textbox_smiles.value = c.isomeric_smiles
            display(textbox_smiles)
            
            textbox_db.value = str(doublebonds)
            display(textbox_db)
            
            textbox_tb.value = str(triplebonds)
            display(textbox_tb)
            
            textbox_rings.value = str(rings)
            display(textbox_rings)

        # Display 2d structure
        with output2:
            try:
                # Get PNG image from PubChem
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{c.cid}/PNG"
                r = requests.get(url)
                r.raise_for_status()
                display(Image(r.content))
            except:
                pass

# Call function when button is clicked
button_getcompound.on_click(get_compound_from_string)

# Set the layout of the output
output_layout = widgets.HBox([output, output2])
display(output_layout)


# In[3]:


# Second button
button_getparameters = widgets.Button(
description='Calculate Parameters',
disabled=False,
)

output3 = widgets.Output()
display(button_getparameters, output3)

# Textboxes for the resulting parameters
textbox_va = widgets.Text(
    description='VA (Schroeder-Volumen):',
    disabled=False   
)

textbox_klipw = widgets.Text(
    description='K lip/w:',
    disabled=False   
)

textbox_dparalip = widgets.Text(
    description='D para lip (cm^^2/s):',
    disabled=False   
)

textbox_dperplip = widgets.Text(
    description='D perp lip (cm^^2/s):',
    disabled=False   
)

textbox_kcorw = widgets.Text(
    description='K cor/w:',
    disabled=False   
)

textbox_dcor = widgets.Text(
    description='D cor (cm^^2/s):',
    disabled=False   
)

# Calculate Parameters
def get_parameters(b):
    output3.clear_output()
    with output3:
        # va
        # Store the atom counts in a dictionary
        atomcounts = {
            "C": 0,
            "H": 0,
            "O": 0,
            "N": 0,
            "Br": 0,
            "Cl": 0,
            "F": 0,
            "I": 0,
            "S": 0
        }
    
        for key in atomcounts:
            atomcounts[key] = textbox_smiles.value.count(key)
    
        if "H" in textbox_formula.value:
            index_h = textbox_formula.value.index("H")
            hcounts = ""
            if not textbox_formula.value[index_h + 1].isdigit():
                hcounts = "1"
            i = 1 
            while textbox_formula.value[index_h + i].isdigit():
                hcounts += textbox_formula.value[index_h + i]
                i += 1
            atomcounts["H"] = int(hcounts)
        #print("Atom counts:\n", atomcounts)
     
        va = 7 * (atomcounts["C"] + atomcounts["H"] + atomcounts["O"] + atomcounts["N"] \
                  + int(textbox_db.value) + 2 * int(textbox_tb.value) - int(textbox_rings.value)) \
            + 31.5 * atomcounts["Br"] + 24.5 * atomcounts["Cl"] + 10.5 * atomcounts["F"] \
            + 38.5 * atomcounts["I"] + 21 * atomcounts["S"]
        
        textbox_va.value = str(va)
        display(textbox_va)
    
        #klipw
        epsilon = 0.00000003
        #print("epsilon =", epsilon)
    
        klipw_intrinsic = 0.43 * (10 ** float(textbox_xlogp.value)) ** 0.81
        #print("klipw_intrinsic =", klipw_intrinsic)
    
        if va <= 342.3:
            a_s = 0.155 * va ** 0.6
        else:
            a_s = 0.735 * va ** (1/3)
        #print("a_s =", a_s)
        
        r_pore = 16
        #print("r_pore =", r_pore)
        
        lambda_microporous = a_s / r_pore
        #print("lambda_microporous =", lambda_microporous)
        
        klipw_microporous = (1 - lambda_microporous) ** 2
        #print("klipw_microporous =", klipw_microporous)
        
        klipw = (1 - epsilon) * klipw_intrinsic + epsilon * klipw_microporous
        #print("klipw =", klipw)
        
        textbox_klipw.value = str(klipw)
        display(textbox_klipw)
     
        # dparalip
        h_para = 1
        #print("h_para =", h_para)
        
        dparalip = (0.000000124 * (100 / float(textbox_weight.value)) ** 2.43 + 0.00000000234) / h_para
        #print("dparalip =", dparalip)
            
        textbox_dparalip.value = str(dparalip)
        display(textbox_dparalip)
    
        # dperplip
        h_perp = 1
        #print("h_perp =", h_perp)
        kperplip_intrinsic = 10 ** (-0.725 - 0.792 * float(textbox_weight.value) ** (1/3)) / h_perp
        #print("kperplip_intrinsic =", kperplip_intrinsic)
        dperplip_intrinsic = kperplip_intrinsic * 0.0000013
        #print("dperplip_intrinsic =", dperplip_intrinsic)
    
        if va <= 342.3:
            d_aq = 0.000188 / va ** 0.6
        else:
            d_aq = 0.0000398 / va ** (1/3)
        #print("d_aq =", d_aq)
        
        dperplip_microporous = d_aq * (1 - 2.104 * lambda_microporous + 2.09 * lambda_microporous ** 3 \
                                       - 0.95 * lambda_microporous ** 5)
        #print("dperplip_microporous =", dperplip_microporous)
        dperplip = ((1 - epsilon) * klipw_intrinsic * dperplip_intrinsic + \
                    epsilon * klipw_microporous * dperplip_microporous) / klipw
        #print("dperplip =", dperplip)
        
        textbox_dperplip.value = str(dperplip)
        display(textbox_dperplip)
    
        #kcorw
        lambda_keratin = a_s / 35
        #print("lambda_keratin =", lambda_keratin)
        phi_f = 0.1928
        #print("phi_f =", phi_f)
        phi_f_strich = phi_f * (1 + lambda_keratin) ** 2
        #print("phi_f_strich =", phi_f_strich)
        kcorw = 1 - phi_f_strich
        #print("kcorw =", kcorw)
        
        textbox_kcorw.value = str(kcorw)
        display(textbox_kcorw)
    
        #dcor
        dcor = d_aq * (1 - phi_f_strich) * \
                (0.9999 - 1.2762 * lambda_keratin + 0.0718 * lambda_keratin ** 2 + 0.1195 * lambda_keratin ** 3)
        #print("dcor =", dcor)
        
        textbox_dcor.value = str(dcor)
        display(textbox_dcor)

# Call function when button is clicked
button_getparameters.on_click(get_parameters)


# In[ ]:





# In[ ]:




