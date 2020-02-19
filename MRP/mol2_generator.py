# -*- coding: utf-8 -*-
"""
This module writes mol2 files with RESP fitted charges for each modified residue. 
"""
from msmtmol.cal import calc_bond
from msmtmol.element import CoRadiiDict

def find_attype(project_name, gatms):
    j=0
    mol2_file = open(project_name + ".mol2", 'r')
    line= mol2_file.readlines()
    while j < len(line):
        if "@<TRIPOS>ATOM" in line[j]:
            start = j +1
        if "@<TRIPOS>BOND" in line[j]:
            end = j
        j +=1
    for k in range (start, end):
        attype = (line[k][50] + line[k][51]).strip()
        gatms[(k-start)][6] = attype
    mol2_file.close()
    return gatms

def parse_for_mol2(gatms, residues, external, id_list):
    molecule_list = residues + external
    parsed_list = []
    parsed_gatms = [0]*len(molecule_list)
    for i in molecule_list:
        parsed_gatms[molecule_list.index(i)] = []    
    k=0
    for i in molecule_list:
        k=0
        l=0
        for j in gatms:
            if j[5] == i:
                parsed_gatms[molecule_list.index(i)].append(j)
        
        parsed_list = parsed_gatms[molecule_list.index(i)]
        bond_list = get_bonded(parsed_list)
        input_file = open(i[:2] + "1.mol2", 'w')
        input_file.write("@<TRIPOS>MOLECULE" + '\n')
        input_file.write(i[:2] + "1" + '\n')
        input_file.write('%5d%6d%6d%6d%6d' %(len(parsed_list), len(bond_list),
                                       1, 0, 0) + '\n')
        input_file.write("SMALL" + '\n')
        input_file.write("RESP Charge" + '\n')
        input_file.write('\n')
        input_file.write('\n')
        input_file.write("@<TRIPOS>ATOM" + '\n')
        while k < len(parsed_list):  
            if parsed_list[k][8] in id_list:
                input_file.write('%7d %-4s    %10.4f%10.4f%10.4f %-4s %6d %-4s %12.6f'\
                        %(k+1, parsed_list[k][4], parsed_list[k][1], parsed_list[k][2], parsed_list[k][3], id_list[parsed_list[k][8]],
                         1, parsed_list[k][5][:2] + "1", float(parsed_list[k][7])) + '\n')
            elif parsed_list[k][4] == "C":
                input_file.write('%7d %-4s    %10.4f%10.4f%10.4f %-4s %6d %-4s %12.6f'\
                        %(k+1, parsed_list[k][4], parsed_list[k][1], parsed_list[k][2], parsed_list[k][3], "c ",
                         1, parsed_list[k][5][:2] + "1", float(parsed_list[k][7])) + '\n')
            else:
                input_file.write('%7d %-4s    %10.4f%10.4f%10.4f %-4s %6d %-4s %12.6f'\
                        %(k+1, parsed_list[k][4], parsed_list[k][1], parsed_list[k][2], parsed_list[k][3], parsed_list[k][6],
                         1, parsed_list[k][5][:2] + "1", float(parsed_list[k][7])) + '\n')
            k +=1
        input_file.write("@<TRIPOS>BOND" + '\n')
        
        while l < len(bond_list):
            input_file.write('%6d%5d%5d%2d' %(l+1, bond_list[l][0],bond_list[l][1], 1) + '\n')
            l +=1
        input_file.write("@<TRIPOS>SUBSTRUCTURE" + '\n')
        input_file.write("     1 " + i[:2] + "1" + "        1 TEMP               " + "0 ****  ****    0 ROOT" + '\n')
        input_file.close()   

def get_bonded(gatms):
    i=0
    j=0
    bond_list = []
    while i < len(gatms):
        crd_1 = [gatms[i][1], gatms[i][2], gatms[i][3]]
        radius_1 = CoRadiiDict[gatms[i][0]]
        while j < len(gatms):
            if j in range(i, len(gatms)):
            
                crd_2 = [gatms[j][1], gatms[j][2], gatms[j][3]]
                radius_2 = CoRadiiDict[gatms[j][0]]
                rad_tot = radius_1 + radius_2 + 0.40
                length = calc_bond(crd_1, crd_2)
                if (length > 0.1) and (length <= rad_tot):
                    bond_list.append([i+1, j+1])
            j +=1
        j = 0
        i +=1
    return bond_list

def parse_adduct(fname, residues, external):
    adduct_file = open(fname + "_capped.pdb", 'r')
    molecule_list = residues + external
    line= adduct_file.readlines()
    for i in molecule_list:
        ind_file = open(i + ".pdb", 'w')
        
        for k in line:
            if i in k:
                ind_file.write(k)
        ind_file.close()
    adduct_file.close()
    