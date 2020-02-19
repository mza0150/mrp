# -*- coding: utf-8 -*-
"""
Title: MRP.py (Modified Residue Parameterizer)
Author: Patrick Sahrmann
Group: Goodwin
University: Auburn

This program is designed to streamline the process of modified residue parameterization.
This includes amino acid residues modified via inhibitor binding, formation of covalent adducts, etc.
MRP.py is designed to work with AmberTools and the Gaussian software package to handle molecular 
modeling and QM calculations, respectively.

"""
from msmtmol.readpdb import get_atominfo_fpdb
from MRP.capped_residues import extract_adduct, write_init_mol2, write_fixed_mol2, write_gaussian
from MRP.resp_fitting import get_isomer_hydrogens, write_resp_input_1, \
write_charges_resp_1, write_resp_input_2, resp_fit, assign_charges
import os
from optparse import OptionParser
from MRP.mol2_generator import parse_adduct, find_attype, parse_for_mol2, get_bonded
from MRP.parameter_search import rewrite_pdb, connect_adduct, write_tleap_input, \
find_missing_params_gaff, find_missing_params_parm10
parser = OptionParser()
parser.add_option("-i", dest="inputfile", type='string',
                  help="Input file name")
parser.add_option("-s", "--step", dest="step", type='string',
                  help="Step number")
(options, args) = parser.parse_args()
project_name = "Project"
pdbname = "PDB"
full_directory = os.getcwd()
directory = full_directory[(full_directory.find("home/")-1):]
residues = []
connect_list = []
id_list = {}
external = []
gatms = []
charge = 0
multiplicity = 1
resp_info = []
reslist = []
extlist = []
gaffatoms = []
corr_isomer = {}
needed_gaff = False
protein_dat = "parm10"
water_model = "TIP3P"
gaff2ff = {"c" : "C",
           "c3" : "CT",
           "hn" : "H",
           "n" : "N",
           "o" : "O",
           "os": "OS", 
           "f" : "F",
           "br" : "Br",
           "cl" : "Cl",
           "i" : "I",
           "ss" : "SH",
           "ca" : "CA",
           "c1" : "CZ",
           "oh" : "OH",
           "s2" : "SH",
           "p5" : "P",
           "s" : "S",
           "ha" : "HA",
           "nb" : "N*"
           }
found_resid = False
found_connect = False
found_gaff = False

inputfile = open(options.inputfile, "r")
for i in inputfile:
    commas = [] 
    gaffatoms = []
    old_connect_list = []
    num_str = ""
    mul_str = ""
    atoms = [0,0]
    j=0
    k=0
    if "pdb_name" in i:
        pdbname = i[8:]
        pdbname = pdbname.replace(" ", "")
        pdbname = pdbname.replace("\n", "")
        pdbname = pdbname.replace(".pdb", "")
    if "project_name" in i:
        project_name = i[12:]
        project_name = project_name.replace(" ", "")
        project_name = project_name.replace("\n", "")
    if "water_model" in i:
        water_model = i[11:]
        water_model = water_model.replace(" ", "")
        water_model = water_model.replace("\n", "")
    if "protein_dat" in i:
        protein_dat = i[11:]
        protein_dat = protein_dat.replace(" ", "")
        protein_dat = protein_dat.replace("\n", "")
    if "charge" in i: 
        for character in i:
            if character.isdigit():
                num_str += character
            if character == "-":
                num_str += character
        charge =  int(num_str)
    if "multiplicity" in i: 
        for character in i:
            if character.isdigit():
                mul_str += character
        multiplicity =  int(mul_str)
    if "residue_ids" in i: 
        found_resid = True
        j=10
        k=0
        commas = []
        for j in range(10, (len(i)-1)):
            if i[j] == ",":
                commas.append(j)
            j +=1 
        j = 10
        if commas != []:
            for k in commas:
                if len(commas) == 1:
                    residues.append(int(i[11:k]))
                    residues.append(int(i[(k+1):]))
                else:
                    if k == min(commas):
                        residues.append(int(i[11:k]))
                    if k == max(commas):
                        residues.append(int(i[(k+1):]))
                    else:
                        residues.append(int(i[(k+1):commas[(commas.index(k)+1)]]))        
                j +=1      
        else:
            residues.append(int(i[11:]))
    if "adduct_connection_atoms" in i: 
        found_connect = True
        j=23
        k=0
        commas = []
        connect = []
        for j in range(23, (len(i)-1)):
            if i[j] == ",":
                commas.append(j)
            j +=1
        if len(commas) == 0:
            connect.append(i[23:])
        for k in commas:
            if len(commas) == 1:
                connect.append(i[23:k])
                connect.append(i[(k+1):])
            else:   
                if k == min(commas):
                    connect.append(i[23:k])
                if k == max(commas):
                    connect.append(i[(k+1):])
                else:
                    connect.append(i[(k+1):commas[(commas.index(k)+1)]])   
        for i in connect:
            i = i.replace("\n", "")
            i = i.replace(" ", "")
            bond = i.index("-")
            atoms[0] = int(i[:bond])
            atoms[1] = int(i[(bond+1):])
            atoms = [atoms[0], atoms[1]]
            old_connect_list.extend(atoms) 
        q = 0
        new_connect_list = [0]*len(old_connect_list)
        connect_list = [0]*int((float(len(old_connect_list))/2))
        while q < (len(old_connect_list)-1):
            new_connect_list[q] = [old_connect_list[q], old_connect_list[q+1]]
            q +=2
        q = 0
        while q < (len(connect_list)):
            connect_list[q] = new_connect_list[2*q]
            q +=1
        
    if "adduct_connection_atom_ids" in i: 
        found_gaff = True
        j=26
        k=0
        commas = []
        while j in range(26, (len(i)-1)):
            if i[j] == ",":
                commas.append(j)
            j +=1 
        for k in commas:
            if len(commas) == 1:
                gaffatoms.append(i[26:k])
                gaffatoms.append(i[(k+1):])
            else:
                if k == min(commas):
                    gaffatoms.append(i[26:k])
                if k == max(commas):
                    gaffatoms.append(i[(k+1):])
                else:
                    gaffatoms.append(i[(k+1):commas[(commas.index(k)+1)]])
        for i in gaffatoms:
            i = i.replace("\n", "")
            i = i.replace(" ", "")
            gid = i.index(":")
            key = int(i[:gid])
            gafftype = i[(gid+1):].replace(" ", "")
            id_list[key] = gafftype 
inputfile.close()

if (options.step == "1"):
    if found_resid == False:
        print("Error: Residue numbers for adduct not specified.")
    if found_connect == True and found_gaff == False:
        print("Error: GAFF type not specified for atoms in adduct_connection_atoms.")
    if found_connect == False and found_gaff == True:
        print("Error: atoms in adduct_connection_atoms not specified.")
    if (found_resid and found_connect and found_gaff):
        print("Writing separate pdb containing capped adduct...")
        reslist = extract_adduct(get_atominfo_fpdb(pdbname + ".pdb")[0], residues, external, project_name, gatms, reslist, extlist)[1]
        print("Creating  mol2 to establish gaff atom types...")
        write_init_mol2(project_name)
        write_fixed_mol2(project_name, gatms, found_gaff, id_list, found_connect, connect_list)
        print("Writing gaussian input files...")
        write_gaussian(project_name, gatms, charge, multiplicity)
        os.remove(project_name + "_init.mol2")
        os.remove("ATOMTYPE.INF")
        os.remove("ANTECHAMBER_BOND_TYPE.AC0")
        os.remove("ANTECHAMBER_BOND_TYPE.AC")
        os.remove("ANTECHAMBER_AC.AC0")
        os.remove("ANTECHAMBER_AC.AC")
    if (found_resid and not found_connect and not found_gaff):
        print("Writing separate pdb containing capped adduct...")
        reslist = extract_adduct(get_atominfo_fpdb(pdbname + ".pdb")[0], residues, external, project_name, gatms, reslist, extlist)[1]
        print("Creating  mol2 to establish gaff atom types...")
        write_init_mol2(project_name)
        write_fixed_mol2(project_name, gatms, found_gaff, id_list, found_connect, connect_list)
        print("Writing gaussian input file...")
        write_gaussian(project_name, gatms, charge, multiplicity)
        os.remove(project_name + "_init.mol2")
        os.remove("ATOMTYPE.INF")
        os.remove("ANTECHAMBER_AC.AC0")
        os.remove("ANTECHAMBER_AC.AC")
   
if (options.step == "2"):
    print("Finding missing intra-adduct parameters...")
    os.system("parmchk2 -i %s.mol2  -f mol2  -o %s.frcmod" %(project_name, project_name))
    reslist = extract_adduct(get_atominfo_fpdb(pdbname + ".pdb")[0], residues, external, project_name, gatms, reslist, extlist)[1]
    bond_full = get_bonded(gatms)
    bond_instruct = connect_adduct(gatms, connect_list, found_connect)
    corr_isomer = get_isomer_hydrogens(bond_full, gatms)
    print("Writing step 1 of 2 stage resp input files...")
    resp_info = write_resp_input_1(project_name, charge, gatms)
    print("Writing initial charge file with provided fixed charges on caps...")
    write_charges_resp_1(project_name, gatms)
    print("Writing step 2 of 2 stage resp input files...")
    write_resp_input_2(project_name, charge, gatms, resp_info, corr_isomer)
    print("Performing resp...")
    resp_fit(project_name)
    find_attype(project_name, gatms)
    assign_charges(project_name, gatms)
    parse_adduct(project_name, reslist, external)
    print("Writing mol2 files for each residue of the adduct...")
    parse_for_mol2(gatms, reslist, extlist, id_list)
    print("Writing modified full protein pdb...")
    rewrite_pdb(pdbname, residues)
    print("Writing initial tleap input file...")
    write_tleap_input(pdbname, residues, reslist, extlist, bond_instruct, 0, external,  water_model, project_name, needed_gaff)
    os.system("tleap -f %s_tleap.in" %(pdbname))
    print("Fetching parameters from parm10.dat...")
    find_missing_params_parm10(pdbname, gaff2ff, protein_dat)
    print("Writing final tleap input file...")
    write_tleap_input(pdbname, residues, reslist, extlist, bond_instruct, 1, external,  water_model, project_name, needed_gaff)
    print("Writing fully parameterized, solvated pdb...")
    os.system("tleap -f %s_tleap.in" %(pdbname))
    if os.path.getsize(pdbname + "_solv.parm7") == 0:
        needed_gaff = True
        print("Possibly missing additional parameters. Fetching parameters from gaff.dat...")
        find_missing_params_gaff(pdbname, gaff2ff)
        print("Writing final tleap input file with additional gaff parameters...")
        write_tleap_input(pdbname, residues, reslist, extlist, bond_instruct, 1)
        os.system("tleap -f %s_tleap.in" %(pdbname))
    print("Removing intermediate files...")
    os.remove(project_name + "_cappedmol.pdb")     
    os.remove(project_name + ".mol2")
    os.remove(project_name + ".pdb")  
    print("MRP.py finished. Ran successfully.")