# -*- coding: utf-8 -*-
"""
This module writes a separate pdb and initial mol2 file containing the modified residues with ACE/NME caps,
then writes the multi-step Gaussian input file including geometry optimization and ESP derivation
for RESP fitting.
"""
from msmtmol.cal import calc_bond
from msmtmol.element import (bdld)
import os

def extract_adduct(mol, residue_ids, external_ids, fname, gatms, reslist, extlist):
    ACE = ['C', 'O', 'CA', 'HA', 'CB', 'N', 'HA2', 'HA3']
    NME = ['N', 'H', 'HN', 'CA', 'HA', 'CB', 'C', 'HA2', 'HA3']
    dict_ace = {"CB":"HA2",
            "N":"HA3"}
    dict_nme = {"CB":"HA2",
            "C":"HA3"}
    adduct_file = open(fname + "_capped.pdb", 'w')
    adduct_file.write('BUILD BY MRP.PY' + '\n')
    #record residue names for each residue in adduct
    for i in residue_ids:
        reslist.append(mol.residues[i].resname)
    #record molecule names for each non-residue in adduct
    for i in external_ids:
        reslist.append(mol.residues[i].resname)  
    #convert necessary atoms in original pdb to obtain ACE/NME caps
    #In ACE:
        #CA -> CH3
        #C, O remain the same
        #HA -> HA1 such that during the second stage of RESP (isomerization) methyl group is recognized
        #CB, N -> HA2, HA3
    #In NME:
        #CA -> CH3
        #HN -> C
        #N, H remain the same
        #HA -> HA1 such that during the second stage of RESP (isomerization) methyl group is recognized
        #CB, C -> HA2, HA3
    for i in residue_ids + external_ids:
        if (i in residue_ids):
            if i != 1:
                for j in mol.residues[i-1].resconter:
                    atom = mol.atoms[j].atname
                    if (atom == 'CA'):
                        cacrd = mol.atoms[j].crd
                for j in mol.residues[i-1].resconter:
                    gtype = "ATOM"
                    atom = mol.atoms[j].atname
                    
                    if (atom in ACE):
                        atomid = mol.atoms[j].atid
                        resid = mol.atoms[j].resid
                        element = mol.atoms[j].element
                        if (atom == 'C') or (atom == 'O'):
                                crdx = mol.atoms[j].crd[0]
                                crdy = mol.atoms[j].crd[1]
                                crdz = mol.atoms[j].crd[2]
                                element = atom 
                        elif (atom == 'CA'):
                                atom = 'CH3'
                                crdx = mol.atoms[j].crd[0]
                                crdy = mol.atoms[j].crd[1]
                                crdz = mol.atoms[j].crd[2]
                                element = 'C'
                        elif (atom == 'HA'):
                                atom = 'HA1'
                                crdx = mol.atoms[j].crd[0]
                                crdy = mol.atoms[j].crd[1]
                                crdz = mol.atoms[j].crd[2]
                        elif (atom in ['CB', 'N']):
                                atom = dict_ace[atom]
                                bvec = calc_bond(cacrd, mol.atoms[j].crd)
                                crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                                crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                                crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                                element = 'H'
                        crdx = round(crdx, 3)
                        crdy = round(crdy, 3)
                        crdz = round(crdz, 3)
                        resname = 'ACE'
                        adduct_file.write("%-6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
                                 %(gtype, atomid, atom, resname, 'A', resid, crdx, crdy, crdz, 1.00, 0.00) + '\n' )
                        gatms.append([element, mol.atoms[j].crd[0], mol.atoms[j].crd[1], mol.atoms[j].crd[2], atom, resname, "", 0, atomid, resid])                      
            for j in mol.residues[i].resconter:
                gtype = "ATOM"
                atm = mol.atoms[j]
                atid = atm.atid
                atomid = mol.atoms[j].atid
                atom = mol.atoms[j].atname
                element = mol.atoms[j].element
                if len(atm.atname) == 3:
                    atname = atm.atname
                else:
                    atname = atm.atname.center(4)
                crd = atm.crd
                resid = atm.resid
                resname = mol.residues[resid].resname
                adduct_file.write("%-6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
                                  %(gtype, atid, atname, resname, 'A', resid, crd[0], crd[1], crd[2], 1.00, 0.00) + '\n' )
                gatms.append([element, mol.atoms[j].crd[0], mol.atoms[j].crd[1], mol.atoms[j].crd[2], atom, resname, "", 0, atomid, resid]) 
                
            for j in mol.residues[i+1].resconter:
                gtype = "ATOM"
                atom = mol.atoms[j].atname
                element = mol.atoms[j].element
                if (atom in NME):
                    if (atom == 'N') or (atom == 'H'):
                        crdx = mol.atoms[j].crd[0]
                        crdy = mol.atoms[j].crd[1]
                        crdz = mol.atoms[j].crd[2]
                        element = atom   
                    elif (atom == 'HN'):
                        atom = 'H'
                        crdx = mol.atoms[j].crd[0]
                        crdy = mol.atoms[j].crd[1]
                        crdz = mol.atoms[j].crd[2]
                        cacrd = [crdx, crdy, crdz]
                        element = 'C' 
                    if (atom == 'CA'):
                        atom = 'CH3'
                        crdx = mol.atoms[j].crd[0]
                        crdy = mol.atoms[j].crd[1]
                        crdz = mol.atoms[j].crd[2]
                        cacrd = [crdx, crdy, crdz]
                        element = 'C'
                    elif (atom == 'HA'):
                        atom = 'HA1'
                        crdx = mol.atoms[j].crd[0]
                        crdy = mol.atoms[j].crd[1]
                        crdz = mol.atoms[j].crd[2]
                        cacrd = [crdx, crdy, crdz]
                    elif (atom in ['CB', 'C']):
                        atom = dict_nme[atom]
                        bvec = calc_bond(cacrd, mol.atoms[j].crd)
                        crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                        crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                        crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                        element = 'H'
                    crdx = round(crdx, 3)
                    crdy = round(crdy, 3)
                    crdz = round(crdz, 3)
                    atomid = mol.atoms[j].atid
                    resid = mol.atoms[j].resid
                    resname = 'NME'
                    adduct_file.write("%-6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
                                %(gtype, atomid, atom, resname, 'A', resid, crdx, crdy, crdz, 1.00, 0.00) + '\n' )
                    gatms.append([element, mol.atoms[j].crd[0], mol.atoms[j].crd[1], mol.atoms[j].crd[2], atom, resname, "", 0, atid, resid])  
        if (i in external_ids):   
            for j in mol.residues[i].resconter:
                    atm = mol.atoms[j]
                    atid = atm.atid
                    atom = mol.atoms[j].atname
                    element = mol.atoms[j].element
                    if len(atm.atname) == 3:
                        atname = atm.atname
                    else:
                        atname = atm.atname.center(4)
                    crd = atm.crd
                    resid = atm.resid
                    resname = mol.residues[resid].resname
                    adduct_file.write("%-6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
                                   %(gtype, atid, atname, resname, 'A', resid, crd[0], crd[1], crd[2], 1.00, 0.00) + '\n' )
                    gatms.append([element, mol.atoms[j].crd[0], mol.atoms[j].crd[1], mol.atoms[j].crd[2], atom, resname, "", 0, atid, resid])  
    adduct_file.write('END')
    adduct_file.close()     
    return gatms, reslist, extlist, mol
    
def write_init_mol2(project_name):
    filename = project_name + "_capped.pdb"
    filename_2 = project_name + "_cappedmol.pdb"
    fileinfo = []
    pdbfile = open(filename, 'r')
    lines = pdbfile.readlines()
    if len(project_name) >= 3:
        i = 1
        while  i < len(lines):
            j = lines[i]
            if j != "END":
                new_line = j[:17] + project_name[0] + project_name[1] + project_name[2] + j[20:]
                fileinfo.append(new_line)
            
            else:
                fileinfo.append(j)
            i +=1
    else:
        i = 0
        while  i < len(pdbfile.readlines()):
            j = lines[i]
            project_name_1 =  project_name + 1 + 1 +1
            if j != "END":
                new_line = j[:17] + project_name_1[0] + project_name_1[1] + project_name_1[2] + j[20:]
                fileinfo.append(new_line)
            else:
                fileinfo.append(j)
            i +=1
    pdbfile.close()
    pdbfile_2 = open(filename_2, 'w')
    for i in fileinfo:
        pdbfile_2.write(i) 
    pdbfile_2.close()
    os.system("antechamber -fi pdb -fo mol2 -i %s -o %s_init.mol2 -j 5 -at gaff -dr no" %(filename_2, project_name))

def write_fixed_mol2(project_name, gatms, found_gaff, id_list, found_connect, connect_list):
    initmol2 = open(project_name + "_init.mol2", 'r')
    mol2 = open(project_name + ".mol2", 'w')
    start = 0 
    end = 0
    j = 0
    atm_1 = 1
    atm_1 = 2
    temp = ""
    start_a = 0
    keys = []
    if found_gaff:
        keys = id_list.keys()
    lines = initmol2.readlines()
    for i in lines:
        if "@<TRIPOS>ATOM" in i:
            start_a = lines.index(i)
        if "@<TRIPOS>BOND" in i:
            start = lines.index(i)
        if "@<TRIPOS>SUBSTRUCTURE" in i:
            end = lines.index(i)
    for k in range(0, (start_a+1)):
        mol2.write(lines[k])
    k = (start_a +1)
    while k in range((start_a+1), start):
        l = k - start_a -1 
        if found_gaff:
            if gatms[l][8] in keys:       
                if len(id_list[gatms[l][8]]) == 1:
                    mol2.write(lines[k][:50] + id_list[gatms[l][8]] + "  " + lines[k][53:])
                else:
                    mol2.write(lines[k][:50] + id_list[gatms[l][8]] + "  " + lines[k][54:])
            else:
                mol2.write(lines[k])
        else:
            mol2.write(lines[k])
        k +=1
    while k in range(start, end):
        mol2.write(lines[k])
        k +=1
    for k in range(end, len(lines)):
        temp += lines[k]
    if found_connect:
        while j < len(connect_list):
            num = (end - start) + 1 + j
            for i in gatms:
                if i[8] == connect_list[j][0]:
                    atm_1 = gatms.index(i)
                if i[8] == connect_list[j][1]:
                    atm_2 = gatms.index(i)
            if num > 10:
                mol2.write("    "  + str(num-1) )
            elif num < 10:
                mol2.write("   "  + str(num-1))
            if atm_1 > 10:
                mol2.write("    " + str(atm_1+1) )
            elif atm_1 < 10:
                mol2.write("   "  + str(atm_1+1) )
            if atm_2 > 10:
                mol2.write("    " + str(atm_2+1) + " 1" + '\n')
            elif atm_2 < 10:
                mol2.write("   "  + str(atm_2+1) + " 1" + '\n')
            j +=1        
    mol2.write(temp)        
    mol2.close()
    
#write multi-step gaussian geometry optimization file and esp 
def write_gaussian(fname, gatms, charge, multip):
    input_file = open(fname + "_esp.com", 'w')
    input_file.write("$RunGauss" + '\n'
               + "%%Chk=%s_opt.chk" %fname +'\n'
               + "%Mem=3000MB" + '\n'
               + "%NProcShared=2" + '\n'
               + "#N B3LYP/6-31G* Geom=PrintInputOrient Integral=(Grid=UltraFine) Opt" + '\n'
               + "SCF = XQC IOp(6/33=2) " + '\n'+ '\n'
               + "Geometry Optimization" + '\n'+ '\n'
               +  "%d %d" %(charge, multip) + '\n')
    i = 0
    while i < len(gatms):
            input_file.write("%-6s %9.4f %9.4f %9.4f" %(gatms[i][0], gatms[i][1], gatms[i][2], gatms[i][3])+ '\n')
            i += 1
    input_file.write('\n')
    input_file.write("--Link1--" + '\n' + '\n'
               + "%%Chk=%s_esp.chk" %fname + '\n'
               + "%%Oldchk=%s_opt.chk" %(fname) + '\n'
               + "%Mem=3000MB" + '\n'
               + "%NProcShared=2" + '\n'
               + "#N B3LYP/6-31G* Geom=check Pop(MK,ReadRadii)" + '\n'
               + "Integral=(Grid=UltraFine) SCF = XQC IOp(6/33=2) Iop(6/42=6)" + '\n'+ '\n'
               + "ESP Derivation" + '\n'+ '\n'
               +  "%d %d" %(charge, multip) + '\n')
    input_file.close
