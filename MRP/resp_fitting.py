# -*- coding: utf-8 -*-
"""
This model generates the RESP input files and performs RESP fitting to obtain charges for the modified 
residues.
"""
from msmtmol.gauio import (get_esp_from_gau)
from periodic_table import AtomicNum
import os

#write RESP input file with restraints on ACE, NME caps
def write_resp_input_1(fname, charge, gatms):
    capnames = ["ACE", "NME"]
    ACE = ['C', 'O', 'CH3', 'HA1', 'HA2', 'HA3']
    NME = ['N', 'H', 'CH3', 'HA1', 'HA2', 'HA3']
    resp_info = []
    i=0
    j=0
    while i < len(gatms):
        element = gatms[i][0]
        atmnum = AtomicNum[element]
        if gatms[i][4] in ACE + NME and gatms[i][5] in capnames:
            resp_info.append([atmnum, -1])
        else:
            resp_info.append([atmnum, 0])
        i +=1
    
    input_file = open(fname + '_resp_charge.in', 'w')
    input_file.write("Resp charges for organic molecule" '\n'
                     + '\n'
                     + " &cntrl" + '\n'
                     + '\n' 
                     + " nmol = 1," + '\n'
                     + " ihfree = 1," + '\n'
                     + " ioutopt = 1," + '\n'
                     + " iqopt = 2," + '\n'
                     + '\n' 
                     + " &end" + '\n'
                     + "    1.0" + '\n'
                     + "Resp charges for organic molecule" + '\n' 
                     + "%5d " %charge + "%4d" %len(gatms) + '\n')
    while j < len(gatms):
        
            if (resp_info[j][0] < 10):
                if (resp_info[j][1] >= 0):
                    input_file.write("    %d    %d" %(resp_info[j][0], resp_info[j][1]) + '\n')
                else:
                    input_file.write("    %d   %d" %(resp_info[j][0], resp_info[j][1]) + '\n')
            else:
                if (resp_info[j][1] >= 0):
                    input_file.write("   %d    %d" %(resp_info[j][0], resp_info[j][1]) + '\n')
                else:
                    input_file.write("   %d   %d" %(resp_info[j][0], resp_info[j][1]) + '\n')
            j += 1
    input_file.write( '\n' )
    input_file.close       
    return resp_info      

#write initial qin file
def write_charges_resp_1(fname, gatms):
    charges_ace = {"C" : 0.597200,
            "O" : -0.567900,
            "CH3" : -0.366200,
            "HA1" : 0.112300,
            "HA2" : 0.112300,
            "HA3" : 0.112300}
    charges_nme = {"N" : -0.415700,
            "H" : 0.271900,
            "CH3" : -0.149000,
            "HA1" : 0.097600,
            "HA2" : 0.097600,
            "HA3" : 0.097600}
    i=0
    null = 0.000000
    input_file = open(fname + '_resp_charge.qin', 'w')
    while i < len(gatms) :
        if ((i+1) % 8 == 0) :
            if (gatms[i][5] == 'ACE'):
                if (gatms[i][4] == 'O') or (gatms[i][4] == 'CH3'):
                    input_file.write(" %f"  %charges_ace[gatms[i][4]] )
                    input_file.write( '\n' )
                else:
                    input_file.write("  %f"  %charges_ace[gatms[i][4]])
                    input_file.write( '\n' )
            elif (gatms[i][5] == 'NME'):
                if (gatms[i][4] == 'N') or (gatms[i][4] == 'CH3'):
                    input_file.write(" %f"  %charges_nme[gatms[i][4]])
                    input_file.write( '\n' )
                else:
                    input_file.write("  %f"  %charges_nme[gatms[i][4]])
                    input_file.write( '\n' )
            else:
                input_file.write("  %f" %null + '\n' )
        else:
            if (gatms[i][5] == "ACE"):     
                if (gatms[i][4] == "O") or (gatms[i][4] == "CH3"):
                    input_file.write(" %f"  %charges_ace[gatms[i][4]])
                else:
                    input_file.write("  %f"  %charges_ace[gatms[i][4]])
            elif (gatms[i][5] == "NME"):
                if (gatms[i][4] == "N") or (gatms[i][4] == "CH3"):
                    input_file.write(" %f"  %charges_nme[gatms[i][4]])
                else:
                    input_file.write("  %f"  %charges_nme[gatms[i][4]])
            else:
                input_file.write("  %f" %null)
        i +=1  
    input_file.close()     

#make the second resp input file in which all methylene/methyl hydrogens are isomerized
def write_resp_input_2(fname, charge, gatms, resp_info, corr_isomer):
    i=0
    l=0
    keys = corr_isomer.keys()
    
    while i < len(gatms):
        if i+1 in keys:
            resp_info[i][1] = (corr_isomer[i+1])
        elif i+1 in corr_isomer.values():
            resp_info[i][1] = 0
        else:
            resp_info[i][1] = -1
        i +=1
    input_file = open(fname + '_resp_h.in', 'w')
    input_file.write("Resp charges for organic molecule" '\n'
                     + '\n'
                     + " &cntrl" + '\n'
                     + '\n' 
                     + " nmol = 1," + '\n'
                     + " ihfree = 1," + '\n'
                     + " ioutopt = 1," + '\n'
                     + " iqopt = 2," + '\n'
                     + '\n' 
                     + " &end" + '\n'
                     + "    1.0" + '\n'
                     + "Resp charges for organic molecule" + '\n' 
                     + "%5d " %charge + "%4d" %len(gatms) + '\n')
    
    while l < len(resp_info):   
            if (resp_info[l][0] < 10):
                if (resp_info[l][1] >= 0):
                    input_file.write("    %d    %d" %(resp_info[l][0], resp_info[l][1]) + '\n')
                else:
                    input_file.write("    %d   %d" %(resp_info[l][0], resp_info[l][1]) + '\n')
            else:
                if (resp_info[l][1] >= 0):
                    input_file.write("   %d    %d" %(resp_info[l][0], resp_info[l][1]) + '\n')
                else:
                    input_file.write("   %d   %d" %(resp_info[l][0], resp_info[l][1]) + '\n')
            l += 1
    input_file.write( '\n' )
    input_file.close              

def assign_charges(project_name, gatms):
    j=0
    mol2_file = open(project_name + "_resp_h.out", 'r')
    line= mol2_file.readlines()
    charges = 0.000000
    while j < len(line):
        if "Point Charges Before & After Optimization" in line[j]:
            start = j +3
        if " Sum over the calculated charges" in line[j]:
            end = j -1
        j +=1
    for k in range (start, end):
        charges = (line[k][30:40]).strip()
        gatms[(k-start)][7] = float(charges)
    mol2_file.close()
    return gatms

def resp_fit(fname):
    resp1 = fname + "_resp_charge"
    resp2 = fname + "_resp_h"
    os.system("espgen -i %s_esp.com.log -o esp.dat" %fname )
    espf = fname + "_esp.dat"
    get_esp_from_gau(fname + "_esp.com.log", espf)
    os.system("resp -O -i %s.in -o %s.out -p %s.pch -q %s.qin -t %s.chq" %(resp1, resp1, resp1, resp1, resp1) \
              + " -e %s -s " %espf + " %s_calc.esp" %resp1)
    os.system("resp -O -i %s.in -o %s.out -p %s.pch -q %s.chq -t %s.chq" %(resp2, resp2, resp2, resp1, resp2) \
              + " -e %s -s " %espf + " %s_calc.esp" %resp2)

def get_isomer_hydrogens(bond_list, gatms):
    isomer_h = []
    isomer_h_copy = []
    h1_list = []
    h2_list = []
    methyl_dict = {}
    methylene_dict = {}
    corr_isomer = {}
    a = 0
    b = 0
    c = 0
    for i in range (0, len(bond_list)):
        for j in range (i+1, len(bond_list)):
            k = bond_list[i]
            l = bond_list[j]
            #if the first entry of both items in the list are equal,
            #that means the second entries of both items are bonded
            #to the same item, therefore check to see if the
            #second entries are hydrogens
            if (k[0] == l[0]):
                #get the element id of the second entry
                h1= gatms[(k[1]-1)][0]
                h2= gatms[(l[1]-1)][0]
                if (h1 == h2) and (h1 == "H"):
                    hydrogens = [k[1], l[1]] 
                    if(k[1] not in isomer_h and l[1] not in isomer_h):
                        isomer_h.append(hydrogens)
            if (k[0] == l[1]):
                h1= gatms[(k[1]-1)][0]
                h2= gatms[(l[0]-1)][0]
                if (h1 == h2) and (h1 == "H"):
                    hydrogens = [k[1], l[0]] 
                    if(k[1] not in isomer_h and l[0] not in isomer_h):
                        isomer_h.append(hydrogens)
            if (k[1] == l[0]):
                h1= gatms[(k[0]-1)][0]
                h2= gatms[(l[1]-1)][0]
                if (h1 == h2) and (h1 == "H"):
                    hydrogens = [k[0], l[1]] 
                    if(k[0] not in isomer_h and l[1] not in isomer_h):
                        isomer_h.append(hydrogens)
            if (k[1] == l[1]):
                h1= gatms[(k[0]-1)][0]
                h2= gatms[(l[0]-1)][0]
                if (h1 == h2) and (h1 == "H"):
                    hydrogens = [k[0], l[0]] 
                    if(k[0] not in isomer_h and l[0] not in isomer_h):
                        isomer_h.append(hydrogens)
    while a <  len(isomer_h):
        h1_list.append(isomer_h[a][0])
        h2_list.append(isomer_h[a][1])
        a +=1
    isomer_h_copy = isomer_h
    k = 0
    while b < len(gatms):
        if (b+1) in h1_list and (b+1) in h2_list:
            index = h1_list.index(b+1) - k
            del isomer_h_copy[index]
            k +=1 
            b += 1
        else:
            b += 1
    isomer_h = isomer_h_copy
    while c < len(isomer_h):
        key = isomer_h[c][0]
        if key in methylene_dict:
           methyl_dict[key] = [methylene_dict[key], isomer_h[c][1]]
           del methylene_dict[key]
        else:
           methylene_dict[key] = isomer_h[c][1]
        c +=1
    for i in isomer_h:
        corr_isomer[i[1]] = i[0]
    return corr_isomer
