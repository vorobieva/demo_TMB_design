import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pyrosetta as py
import glob
import os

py.init("-mute all")

# First step: check that the hydrogen bond of the motif is still present after relaxation. We start with 2493 designs.

hbond_set = py.rosetta.core.scoring.hbonds.HBondSet()
pdbs_list = []

for pdb_i in glob.glob("*_input_*.pdb"):
    hbonds = {'acc':[], 'don':[]}
    folder = os.path.splitext(pdb_i)[0]
    with open(pdb_i, 'r') as fi:
        for line in fi:
            if "HBNet" in line:
                if int(line.split()[2]) in [26,114] :
                    hbonds['don'].append(int(line.split()[2]))
                else:
                    hbonds['acc'].append(int(line.split()[2]))
    t = len(hbonds['don'])*10
    n = 0
    for i in range(1,11):
        design = folder + '/' + folder + '_' + "{:04d}".format(i) + '.pdb'
        pose = py.pose_from_pdb(design)
        pose.update_residue_neighbors()
        py.rosetta.core.scoring.hbonds.fill_hbond_set(pose, False, hbond_set)
        for hb in range(1, hbond_set.nhbonds()+1):
            if hbond_set.hbond(hb).don_res() in hbonds['don'] and hbond_set.hbond(hb).acc_res() in hbonds['acc'] and hbond_set.hbond(hb).don_hatm_is_backbone() == False and hbond_set.hbond(hb).acc_atm_is_backbone() == False:
                n += 1
    if n/t >= 0.7:
        pdbs_list.append(pdb_i)
#        print(pdb_i)
#print(len(pdbs_list))
with open("pdbs_stable_motifs.csv", 'w') as out_f:
    for i in pdbs_list:
        out_f.write(i+"\n")
