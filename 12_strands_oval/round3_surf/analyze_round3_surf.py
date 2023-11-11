import pandas as pd
import numpy as np
import glob
import os
import re
from shutil import copyfile
from collections import Counter
import pyrosetta as py
py.init()

def expected_hbonds(df):
	hbond_dict = {}
	for index,row in df.iterrows():
		hbonds = []
		cst = row["network_cst"]
		lines = cst.split('\n')
		for line in lines:      
			if line != "" and "AtomPair " in line:
				tyr_10 = False
				tyr_96 = False            
				acc_res = int(line.split()[2])
				don_res = int(line.split()[4])
				acc_atm = line.split()[1]
				don_atm = line.split()[3]
				if " OH  10 " in line:
					tyr_10 = True
				if " OH  96 " in line:
					tyr_96 = True         
				if acc_atm not in ["CB","N","CA","O"] and don_atm not in ["CB","N","CA","O"]:
					hbonds.append((acc_res,don_res,acc_atm,don_atm,tyr_10,tyr_96))
		hbond_dict[index] = hbonds
	return (hbond_dict)

def calculate_retention(df,hbond_dict):
	hbond_set = py.rosetta.core.scoring.hbonds.HBondSet()
	dretained = pd.DataFrame(columns=['%_retained','%_retained_tyr'])

	for index,row in df.iterrows():
		interactions = hbond_dict[index]
		interact_res = Counter([inter[0:2] for inter in interactions])
		interact_res_tyr = Counter([inter[0:2] for inter in interactions if True in inter])
		num_interactions = float(len(list(interact_res.elements())))
		num_interactions_tyr = float(len(list(interact_res_tyr.elements())))
		curr_retention = []
		curr_retention_tyr = []
		for pdb in glob.glob(index + "/" + index + "_*.pdb"):
			current_inter_list = []
			pose = py.pose_from_pdb(pdb)
			pose.update_residue_neighbors()
			py.rosetta.core.scoring.hbonds.fill_hbond_set(pose, False, hbond_set)
			for i in range(1, hbond_set.nhbonds()+1):
				if hbond_set.hbond(i).acc_atm_is_backbone() == False and hbond_set.hbond(i).don_hatm_is_backbone() == False:
					current_inter_list.append((hbond_set.hbond(i).acc_res(),hbond_set.hbond(i).don_res()))
					current_inter = Counter(current_inter_list)
					num_curr_inter = 0.0
					num_curr_inter_tyr = 0.0
					for inter in current_inter:
						if inter in interact_res:
							if current_inter[inter] >= interact_res[inter]:
								num_curr_inter += interact_res[inter]
							elif current_inter[inter] < interact_res[inter]:
								num_curr_inter += current_inter[inter]

						if inter in interact_res_tyr:
							if current_inter[inter] >= interact_res_tyr[inter]:
								num_curr_inter_tyr += interact_res_tyr[inter]
							elif current_inter[inter] < interact_res_tyr[inter]:
								num_curr_inter_tyr += current_inter[inter]
					curr_retention.append(num_curr_inter/num_interactions)
					curr_retention_tyr.append(num_curr_inter_tyr/num_interactions_tyr)
		dretained.loc[index] = [sum(curr_retention)/len(curr_retention), sum(curr_retention_tyr)/len(curr_retention_tyr)]
	return(dretained)

if __name__ == "__main__":
	df = pd.read_pickle('../picked_networks.pickle')
	hbond_dict = expected_hbonds(df)
	dretained = calculate_retention(df,hbond_dict)
	dretained.to_pickle('networks_hbond_retention.pickle')
