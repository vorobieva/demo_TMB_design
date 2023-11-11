import pickle

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pyrosetta as py

import glob
import os

py.init("-mute all")

pdbs_list = []

with open("pdbs_stable_motifs.csv", 'r') as in_f:
    for line in in_f:
        pdbs_list.append(line.strip())

energy_diff = {}
for pdb in pdbs_list:
    folder = os.path.splitext(pdb)[0]
    HBNet = []
    energy_table_A = []
    tot_energy_A = 0
    mean_energy_A = 0
    tot_hbond_A = 0
    mean_hbond_A = 0
    labels = []
    with open(pdb, 'r') as pdb_f:
        for line in pdb_f:
            if "HBNet" in line and int(line.split()[2]) != 26 and int(line.split()[2]) != 114 :
                HBNet.append(int(line.split()[2]))
            elif line.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
                labels = pdb_f.readline().strip().split(' ')
                next(pdb_f)
                next(pdb_f)
                for line in pdb_f:
                    if line.startswith("#END_POSE_ENERGIES_TABLE"):
                        break
                    energy_table_A.append(line.strip().split(' '))
    da = pd.DataFrame(energy_table_A, columns=labels)
    da['label'] = da['label'].str.replace('HIS_D', 'HIS')
    new_c = da['label'].str.split('_', expand=True)
    da['resn'] = new_c[1].astype(int) 
    tot_energy_A = da['total'].loc[da['resn'].isin(HBNet)].astype(float).sum()
    mean_energy_A = tot_energy_A/len(HBNet)
    tot_hbond_A = da['hbond_sc'].loc[da['resn'].isin(HBNet)].astype(float).sum()
    mean_hbond_A = tot_hbond_A/len(HBNet)    
    
    for i in range(1,11):
        design = folder + '/' + folder + '_' + "{:04d}".format(i) + '.pdb'
        tot_energy_B = 0
        mean_energy_B = 0
        tot_hbond_B = 0
        mean_hbond_B = 0
        energy_table_B = []    
        diff_energy = 0
        diff_hbond = 0
        with open(design, 'r') as fi:
            for line in fi:
                if line.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
                    next(fi)
                    next(fi)
                    next(fi)
                    for line in fi:
                        if line.startswith("#END_POSE_ENERGIES_TABLE"):
                            break
                        energy_table_B.append(line.strip().split(' '))
        db = pd.DataFrame(energy_table_B, columns=labels)
        db['label'] = db['label'].str.replace('HIS_D', 'HIS')
        new_c = db['label'].str.split('_', expand=True)
        db['resn'] = new_c[1].astype(int)
        tot_energy_B = db['total'].loc[db['resn'].isin(HBNet)].astype(float).sum()
        mean_energy_B = tot_energy_B/len(HBNet)
        tot_hbond_B = db['hbond_sc'].loc[db['resn'].isin(HBNet)].astype(float).sum()
        mean_hbond_B = tot_hbond_B/len(HBNet)
        diff_energy = mean_energy_A - mean_energy_B
        diff_hbond = mean_hbond_A - mean_hbond_B
        energy_diff[design] = (diff_energy, diff_hbond)
de = pd.DataFrame.from_dict(energy_diff, orient='index', columns=['diff_energy','diff_hbond'])

de['design'], de['surface'] = de.index.str.split('/',1).str

dmean = de.groupby('design').mean()

de["diff_sum"] = de['diff_energy'] + de['diff_hbond']
de["rank"] = de.groupby('design')["diff_sum"].rank(ascending=True)

hbond_set = py.rosetta.core.scoring.hbonds.HBondSet()
designs_list = []

for index, row in dmean.iterrows():
    model = index
#    parent = '_'.join(index.split('_')[:-1])
#    print(parent)
#    score_file = model + "/" + "score.sc"
#    ds = pd.read_csv(score_file, sep='\s+', header=1)
#    cst_i = '../' + parent + '/cst'
    hbonds = {'acc':[], 'don':[]}

    good_designs = []
    with open(model+".pdb", 'r') as fi:
        for line in fi:
            if "HBNet" in line:
                if int(line.split()[2]) in [26,114] :
                    hbonds['don'].append(int(line.split()[2]))
                else:
                    hbonds['acc'].append(int(line.split()[2]))
    t = len(hbonds['don'])

    for i in range(1,11):
        n = 0
        design = model + '/' + model + '_' + "{:04d}".format(i) + '.pdb'
        pose = py.pose_from_pdb(design)
        pose.update_residue_neighbors()
        py.rosetta.core.scoring.hbonds.fill_hbond_set(pose, False, hbond_set)
        for hb in range(1, hbond_set.nhbonds()+1):
            if hbond_set.hbond(hb).don_res() in hbonds['don'] and hbond_set.hbond(hb).acc_res() in hbonds['acc'] and hbond_set.hbond(hb).don_hatm_is_backbone() == False and hbond_set.hbond(hb).acc_atm_is_backbone() == False:
                n += 1
        if n == t :
            good_designs.append(model+'_'+"{:04d}".format(i)+'.pdb')
    print(good_designs)

    top_found = False
    r = 1
    while top_found == False and r < 11:
        print(r)
        top = de.loc[(de["rank"] == r) & (de["design"] == model)].surface.item()
        if top in good_designs:
            top_found = True
            designs_list.append(model+"/"+top)
        else:
            r += 1


with open('final_picked_pdbs.pkl', 'wb') as f:
    pickle.dump(designs_list, f)



