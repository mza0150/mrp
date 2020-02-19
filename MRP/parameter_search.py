# -*- coding: utf-8 -*-
"""
This module writes the necessary frcmod files to fully parameterize the adduct, 
as well as the tleap input file. Once fully parameterized, the module runs tleap to obtain a 
solvated protein ready for minimization.
"""

#write commands necessary to connect molecules within adduct
def connect_adduct(gatms, connect_list, found_connect):
    bond_instruct = ""
    connect_dict = {}
    for i in gatms:
        connect_dict[i[8]] = [i[9], i[4]]
    if found_connect:
        for j in connect_list:
            bond_instruct += "bond prot." + str(connect_dict[[j][0][0]][0]) + "." + str(connect_dict[[j][0][0]][1])
            bond_instruct += " prot." + str(connect_dict[[j][0][1]][0])+ "." + str(connect_dict[[j][0][1]][1])  + '\n'   
    return bond_instruct

def rewrite_pdb(pdbname, residues):
    pdb_init = open(pdbname + ".pdb", 'r')
    pdb_final = open(pdbname + "_MRP.pdb", 'w')
    for i in pdb_init.readlines():
        res = i[23:26]
        resid = res.strip()
        if "END" not in i and "CONECT" not in i:
            if "MODEL" not in i:
                if int(resid) in residues:
                    pdb_final.write(i[:19] + "1" + i[20:])
                else:
                    pdb_final.write(i)
        else:
            pdb_final.write("END")
    pdb_init.close()
    pdb_final.close()

def write_tleap_input(pdbname, residues, reslist, extlist, bond_instruct, missingparam, external, water_model, project_name, needed_gaff):
    molecule_list = reslist + external
    tleap = open(pdbname + "_tleap.in", "w")
    tleap.write("source leaprc.protein.ff14SB" + '\n')
    tleap.write("source leaprc.gaff" + '\n')
    tleap.write("source leaprc.water." + water_model.lower() + '\n')
    tleap.write("loadamberparams " + project_name + ".frcmod" + '\n')
    if missingparam == 1:
        tleap.write("loadamberparams " + pdbname + "_MRP.frcmod" + '\n')
    if needed_gaff:
        tleap.write("loadamberparams " + pdbname + "_gaff_MRP.frcmod" + '\n')
    for i in molecule_list:
        tleap.write(i[:2] + "1 = loadmol2 " + i[:2] + "1.mol2" + '\n')
    tleap.write("prot = loadpdb " + pdbname + "_MRP.pdb" + '\n')
    for i in residues:
        if i !=1:
            tleap.write("bond prot." + str(i-1) + ".C " + "prot." + str(i) + ".N" + '\n')
        tleap.write("bond prot." + str(i) + ".C " + "prot." + str(i+1) + ".N" + '\n')
    tleap.write(bond_instruct)
    tleap.write("savepdb prot " + pdbname + "_dry.pdb "  + '\n')
    tleap.write("saveamberparm prot " + pdbname + "_dry.parm7 " +  pdbname + "_dry.rst7 " + '\n')
    if missingparam == 1:
        tleap.write("solvateoct prot " + water_model.upper() + "BOX 20.0" + '\n')
        tleap.write("addions prot Na+ 0" + '\n')
        tleap.write("addions prot Cl- 0" + '\n')
        tleap.write("savepdb prot " +  pdbname + "_solv.pdb"  + '\n')
        tleap.write("saveamberparm prot " + pdbname + "_solv.parm7 " +  pdbname + "_solv.rst7 " + '\n')
    tleap.write("quit")
    tleap.close()
    
def find_missing_params_parm10(pdbname, gaff2ff, protein_dat):
    tleap = open("leap.log", "r")
    ff = open(protein_dat + ".dat", "r")
    frcmod = open(pdbname + "_MRP.frcmod", "w")
    angle_start = 0
    torsion_start = 0
    im_torsion_start = 0
    i = 0
    j = 0
    k = 0
    l = 0
    m = 0
    partition = []
    ff_angle_list = []
    ff_torsion_list = []
    angle_params = []
    torsion_params = []
    write_angle_params = []
    write_torsion_params = []
    log_text = tleap.readlines()
    parm99_text = ff.readlines()
    while i < len(log_text):
        if "Building angle parameters" in log_text[i]:
            angle_start = i
        if "Building proper torsion parameters." in log_text[i]:
            torsion_start = i
        if "Building improper torsion parameters." in log_text[i]:
            im_torsion_start = i
        i +=1
        
    while m < len(parm99_text):   
        if parm99_text[m] == "\n":
            partition.append(m)
        m +=1     
    for k in range (angle_start+1, torsion_start):
        new_angle_text = ""
        original_angle_text = ""
        l = 0
        m = 0
        angle_text = log_text[k]
        angle_text = angle_text[32:]
        if "Error" not in angle_text and "" != angle_text:
            angle_text = angle_text.replace('\n', "")           
            gaffats = [[angle_text[0], 1], [angle_text[4], 5], [angle_text[8], 9]]
            while m < len(gaffats):
                i = gaffats[m]
                if i[1] == len(angle_text):
                    gaffat = i[0] 
                else:
                    gaffat = i[0] + angle_text[i[1]]
                gaffat = gaffat.replace(" ", "")
                if len(gaffat) == 2:
                    original_angle_text += gaffat + "-"
                else:
                    original_angle_text += gaffat + " -"
                original_angle_text = original_angle_text[:8]
                m +=1
            
            while l < len(gaffats):
                i = gaffats[l]
                if i[1] == len(angle_text):
                    gaffat = i[0] 
                else:
                    gaffat = i[0] + angle_text[i[1]]
                gaffat = gaffat.replace(" ", "")
                if gaffat[0].islower():
                    gaffat = gaff2ff[gaffat]
                if len(gaffat) == 2:
                    new_angle_text += gaffat + "-"
                else:
                    new_angle_text += gaffat + " -"
                    new_angle_text = new_angle_text[:8]
                l +=1
            ff_angle_list.append([original_angle_text, new_angle_text])
        
    for k in range (torsion_start+1, im_torsion_start):
        new_torsion_text = ""
        original_torsion_text = ""
        l = 0
        m = 0
        torsion_text = log_text[k]
        if "Error" not in torsion_text and "\n" != torsion_text:
            torsion_text = torsion_text[26:]
            torsion_text = torsion_text.replace('\n', "")
            while m < len(torsion_text):
                if m == len(torsion_text) -1:
                    original_torsion_text += torsion_text[m] 
                else:
                    if torsion_text[m] != "-" and torsion_text[m+1] != "-":
                        original_torsion_text += torsion_text[m] + torsion_text[m+1]  + "-"
                        m+=1
                    elif torsion_text[m] != "-" and torsion_text[m+1] == "-":
                        original_torsion_text += torsion_text[m] + " -"
                m+=1
            if original_torsion_text[len(original_torsion_text)-1] == "-":
                original_torsion_text = original_torsion_text[:len(original_torsion_text)-1]
        
            while l < len(torsion_text):
                if l == len(torsion_text) -1:
                    new_torsion_text += torsion_text[l] 
                else:
                    if torsion_text[l] != "-" and torsion_text[l+1] != "-":
                        if torsion_text[l].islower():
                            new_torsion_text += gaff2ff[torsion_text[l] + torsion_text[l+1]]  + "-"
                        else: 
                            new_torsion_text += torsion_text[l] + torsion_text[l+1]  + "-"
                        l+=1
                    elif torsion_text[l] != "-" and torsion_text[l+1] == "-":
                        if torsion_text[l].islower():
                            new_torsion_text += gaff2ff[torsion_text[l]] + " -"
                        else:
                            new_torsion_text += torsion_text[l] + " -"
                l+=1
            if new_torsion_text[(len(new_torsion_text)-1)] == "-":
                new_torsion_text = new_torsion_text[:len(new_torsion_text)-1]
            ff_torsion_list.append([original_torsion_text, new_torsion_text])
                    
    for i in ff_angle_list:
        if i[1][(len(i[1])-1)] == "-":
            i[1] = i[1][:(len(i[1])-1)]
        i_reverse = i[1][6:] + i[1][2:6] + i[1][:2] 
        for k in range(partition[1],partition[2]):
            if i[1] in parm99_text[k]  or i_reverse in  parm99_text[k]:
                angle_params.append(parm99_text[k][8:])
                write_angle_params.append(i[0] + parm99_text[k][(len(i[1])+1):])
    
    for i in ff_torsion_list:
        i_reverse = ""
        i_X = ""
        i_X_reverse = ""
        j = len(i[1])
        i_reverse += i[1][(j-2):]  + " " + i[1][5:8] + i[1][2:5] + "-"  + i[1][:2]
        i_reverse = i_reverse.lstrip()
        i_X += "X "  + i[1][2:5] + i[1][5:8] + "-X"
        i_X_reverse += "X " + i[1][5:8] + i[1][2:5] + "-X"
        for k in range(partition[2],len(parm99_text)):
            if i[1] in parm99_text[k]  or i_reverse in  parm99_text[k]:
                torsion_params.append(parm99_text[k][11:])
                write_torsion_params.append(i[0] + parm99_text[k][(len(i[1])+1):])
            elif i_X in parm99_text[k]  or i_X_reverse in  parm99_text[k]:
                torsion_params.append(parm99_text[k][11:])
                write_torsion_params.append(i[0] + parm99_text[k][(len(i[1])+1):])

    frcmod.write("PARAMETERS FOUND BY MRP.py"  + '\n')
    frcmod.write("BOND" + '\n')
    frcmod.write("C - n   490.0    1.335       JCC,7,(1986),230; AA" + '\n')
    frcmod.write("c - N   490.0    1.335       JCC,7,(1986),230; AA" + '\n')
    frcmod.write("C - n   490.0    1.335       JCC,7,(1986),230; AA" + '\n')
    frcmod.write('\n' + "ANGLE" + '\n')
    for i in write_angle_params:
        frcmod.write(i)
    frcmod.write('\n' + "DIHE" + '\n')
    for i in write_torsion_params:
        frcmod.write(i)
    frcmod.write('\n' + "IMPR" + '\n' + '\n')
    frcmod.close()
    ff.close()
    tleap.close()

def find_missing_params_gaff(pdbname, gaff2ff):
    tleap = open("leap.log", "r")
    ff = open("gaff.dat", "r")
    frcmod = open(pdbname + "_gaff_MRP.frcmod", "w")
    angle_start = 0
    torsion_start = 0
    im_torsion_start = 0
    i = 0
    j = 0
    k = 0
    l = 0
    m = 0
    partition = []
    ff_angle_list = []
    ff_torsion_list = []
    angle_params = []
    torsion_params = []
    write_angle_params = []
    write_torsion_params = []
    log_text = tleap.readlines()
    parm99_text = ff.readlines()
    while i < len(log_text):
        if "Building angle parameters" in log_text[i]:
            angle_start = i
        if "Building proper torsion parameters." in log_text[i]:
            torsion_start = i
        if "Building improper torsion parameters." in log_text[i]:
            im_torsion_start = i
        i +=1
        
    while m < len(parm99_text):   
        if parm99_text[m] == "\n":
            partition.append(m)
        m +=1     
    for k in range (angle_start+1, torsion_start):
        new_angle_text = ""
        original_angle_text = ""
        l = 0
        m = 0
        angle_text = log_text[k]
        angle_text = angle_text[32:]
        if "Error" not in angle_text and "" != angle_text:
            angle_text = angle_text.replace('\n', "")           
            gaffats = [[angle_text[0], 1], [angle_text[4], 5], [angle_text[8], 9]]
            while m < len(gaffats):
                i = gaffats[m]
                if i[1] == len(angle_text):
                    gaffat = i[0] 
                else:
                    gaffat = i[0] + angle_text[i[1]]
                gaffat = gaffat.replace(" ", "")
                if len(gaffat) == 2:
                    original_angle_text += gaffat + "-"
                else:
                    original_angle_text += gaffat + " -"
                original_angle_text = original_angle_text[:8]
                m +=1
            
            while l < len(gaffats):
                i = gaffats[l]
                if i[1] == len(angle_text):
                    gaffat = i[0] 
                else:
                    gaffat = i[0] + angle_text[i[1]]
                gaffat = gaffat.replace(" ", "")
                if gaffat[0].isupper():
                    gaffat = [key for key, value in gaff2ff.items() if value == gaffat][0]
                if len(gaffat) == 2:
                    new_angle_text += gaffat + "-"
                else:
                    new_angle_text += gaffat + " -"
                    new_angle_text = new_angle_text[:8]
                l +=1
            ff_angle_list.append([original_angle_text, new_angle_text])
        
    for k in range (torsion_start+1, im_torsion_start):
        new_torsion_text = ""
        original_torsion_text = ""
        l = 0
        m = 0
        torsion_text = log_text[k]
        if "Error" not in torsion_text and "\n" != torsion_text:
            torsion_text = torsion_text[26:]
            torsion_text = torsion_text.replace('\n', "")
            while m < len(torsion_text):
                if m == len(torsion_text) -1:
                    original_torsion_text += torsion_text[m] 
                else:
                    if torsion_text[m] != "-" and torsion_text[m+1] != "-":
                        original_torsion_text += torsion_text[m] + torsion_text[m+1]  + "-"
                        m+=1
                    elif torsion_text[m] != "-" and torsion_text[m+1] == "-":
                        original_torsion_text += torsion_text[m] + " -"
                m+=1
            if original_torsion_text[len(original_torsion_text)-1] == "-":
                original_torsion_text = original_torsion_text[:len(original_torsion_text)-1]
        
            while l < len(torsion_text):
                if l == len(torsion_text) -1:
                    new_torsion_text += torsion_text[l] 
                else:
                    if torsion_text[l] != "-" and torsion_text[l+1] != "-":
                        if torsion_text[l].isupper():
                            new_torsion_text += [key for key, value in gaff2ff.items() if value == torsion_text[l]][0] + torsion_text[l+1]  + "-"
                        else: 
                            new_torsion_text += torsion_text[l] + torsion_text[l+1]  + "-"
                        l+=1
                    elif torsion_text[l] != "-" and torsion_text[l+1] == "-":
                        if torsion_text[l].isupper():
                            new_torsion_text += [key for key, value in gaff2ff.items() if value == torsion_text[l]][0] + " -"
                        else:
                            new_torsion_text += torsion_text[l] + " -"
                l+=1
            if new_torsion_text[(len(new_torsion_text)-1)] == "-":
                new_torsion_text = new_torsion_text[:len(new_torsion_text)-1]
            ff_torsion_list.append([original_torsion_text, new_torsion_text])
                    
    for i in ff_angle_list:
        if i[1][(len(i[1])-1)] == "-":
            i[1] = i[1][:(len(i[1])-1)]
        i_reverse = i[1][6:] + i[1][2:6] + i[1][:2] 
        for k in range(partition[1],partition[2]):
            if i[1] in parm99_text[k]  or i_reverse in  parm99_text[k]:
                angle_params.append(parm99_text[k][8:])
                write_angle_params.append(i[0] + parm99_text[k][(len(i[1])+1):])
    
    for i in ff_torsion_list:
        i_reverse = ""
        i_X = ""
        i_X_reverse = ""
        j = len(i[1])
        i_reverse += i[1][(j-2):]  + " " + i[1][5:8] + i[1][2:5] + "-"  + i[1][:2]
        i_reverse = i_reverse.lstrip()
        i_X += "X "  + i[1][2:5] + i[1][5:8] + "-X"
        i_X_reverse += "X " + i[1][5:8] + i[1][2:5] + "-X"
        for k in range(partition[2],len(parm99_text)):
            if i[1] in parm99_text[k]  or i_reverse in  parm99_text[k]:
                torsion_params.append(parm99_text[k][11:])
                write_torsion_params.append(i[0] + parm99_text[k][(len(i[1])+1):])
            elif i_X in parm99_text[k]  or i_X_reverse in  parm99_text[k]:
                torsion_params.append(parm99_text[k][11:])
                write_torsion_params.append(i[0] + parm99_text[k][(len(i[1])+1):])

    frcmod.write("PARAMETERS FOUND BY MRP.py"  + '\n')
    frcmod.write("BOND" + '\n')
    frcmod.write("C - n   490.0    1.335       JCC,7,(1986),230; AA" + '\n')
    frcmod.write("c - N   490.0    1.335       JCC,7,(1986),230; AA" + '\n')
    frcmod.write("C - n   490.0    1.335       JCC,7,(1986),230; AA" + '\n')
    frcmod.write('\n' + "ANGLE" + '\n')
    for i in write_angle_params:
        frcmod.write(i)
    frcmod.write('\n' + "DIHE" + '\n')
    for i in write_torsion_params:
        frcmod.write(i)
    frcmod.write('\n' + "IMPR" + '\n' + '\n')
    frcmod.close()
    ff.close()
    tleap.close()
