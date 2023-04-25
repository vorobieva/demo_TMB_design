import itertools
import re,sys,os,copy
from argparse import ArgumentParser
import numpy as np
from BuildBlueprints import *
from Blueprint import Blueprint

#==============================
# INPUT PARAMETERS
#==============================
parser = ArgumentParser()
parser.add_argument('-xml', type=str, help="xml template")
parser.add_argument('-refpdb', type=str, help="input pdb")
parser.add_argument('-prefix', type=str)
args = parser.parse_args()

template_xml = args.xml

# if there is a bulge in E2, then it should be 13 residues long. If non-bulged, then 12
topol = "E[23]L[3]E[24]L[2]E[22]L[3]E[24]L[2]E[22]L[3]E[26]L[2]E[24]L[3]E[26]L[2]E[22]L[3]E[24]L[2]E[22]L[3]E[24]L[1]"
common_bulges = {'E2': [24], 'E4': [24], 'E6': [26], 'E8': [26], 'E10': [24], 'E12': [24]}
cap = 2
# E2:8,10
# bulges = {}
cwd = os.getcwd()
refpdb = os.path.join(cwd, args.refpdb)
topology = [('E1','E2'),('E2','E3'),('E3','E4'),('E4','E5'),('E5','E6'),('E6','E7'),('E7','E8'),('E8','E9'),('E9','E10'),('E10','E11'),('E11','E12'),('E12','E1')]

ss,combinations,bulges = getCombinationsForBetaBarrels(topol,common_bulges)

# Make directory and bluprints for each combination along with bulge positions
for bulge, comb in zip(bulges, combinations):
#for comb in combinations:
	# Build directories and add bulge combinations
	#---------------------------------------------------------
	# Make the directory name
	comb_name = ""
	for i,s in enumerate(ss):
		comb_name+='%s%i' %(ss[i],comb[i])
	ss_with_corner, comb_with_corner = add_bserine_corner(ss, comb)
	blueprint_name = comb_name+'_bp'	
	cst_name = comb_name+'_cst'

	# First build a plain blueprint file to after creat a bp object to make easier the bulge combinations
#	MakePlainBlueprint(ss,comb,'bp')
	MakeRefBlueprint(ss_with_corner,comb_with_corner,refblue=blueprint_name, bulges=bulge, mostCommonLoopABEGO=True, strandsABEGO="B", cap=cap)
	blue = Blueprint(blueprint_name)
	blue.reindex_blueprint(start=1)

	# copy files provided by the first step
	## Build blueprints
#	MakeRefBlueprint(ss,comb,bulges=bulge,refblue = 'bp')

	#---------------------------------------------------------
	# Blueprints are already written. Now XMLs, cst...

	# names of blueprints must be consistent between here and the template.xml
	# Move above dir for another topology
	################################################
	# step1
	cst_fileout = open(cst_name,'w')
	for hairpin in topology:
#		if hairpin[0] == 'E1':
#			if  cap == 1:
#				sthb = HbondsRegularHairpin(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue,str1_offset=1); cst_fileout.write(sthb)
#			elif cap == 2:
				
		if abs(int(hairpin[1][1:])-int(hairpin[0][1:])) == 1:
			# bulges can not happen between non-sequential strands for now
			n_bulges = 0
			bulge_to_add = None
			if hairpin[0] in bulge:
				# for now, it is only possible to add 1 bulge per hairpin, which should be sufficient for now if we consider the G1 bulge as part of the beta-turn
				bulge[hairpin[0]].sort(reverse=True)
				for b in bulge[hairpin[0]]:
					bulgepos = blue.segment_dict[hairpin[0]].bp_data[b-1][0]
					if (blue.segment_dict[hairpin[0]].bp_data[-1][0] - bulgepos - n_bulges) %2 != 0:
						n_bulges += 1
						bulge_to_add = bulgepos
			elif hairpin[1] in bulge:
				n_bulges = 0
				bulge[hairpin[1]].sort()
				for b in bulge[hairpin[1]]:
					bulgepos = blue.segment_dict[hairpin[1]].bp_data[b-1][0]
					if (bulgepos - blue.segment_dict[hairpin[1]].bp_data[0][0] - n_bulges) %2 == 0:
						n_bulges += 1
						bulge_to_add = bulgepos
			if n_bulges == 1:
				sthb = HbondsBulgedStrand(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue,bulge_position=bulge_to_add) ; cst_fileout.write(sthb)
				cst_fileout.write('\n')
			elif n_bulges == 0:
				sthb = HbondsRegularHairpin(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue) ; cst_fileout.write(sthb)
				cst_fileout.write('\n')
			elif n_bulges > 1:
				print("Something wrong with bulge count!")
		else:
			n_bulges = 0
			bulge_to_add = None
			if hairpin[0] in bulge:
				# for now, it is only possible to add 1 bulge per hairpin, which should be sufficient for now if we consider the G1 bulge as part of the beta-turn
				bulge[hairpin[0]].sort(reverse=True)
				for b in bulge[hairpin[0]]:
					bulgepos = blue.segment_dict[hairpin[0]].bp_data[b-1][0]
					if (blue.segment_dict[hairpin[0]].bp_data[-1][0] - bulgepos - n_bulges) %2 != 0:
						n_bulges += 1
						bulge_to_add = bulgepos
			elif hairpin[1] in bulge:
				n_bulges = 0
				bulge[hairpin[1]].sort()
				for b in bulge[hairpin[1]]:
					bulgepos = blue.segment_dict[hairpin[1]].bp_data[b-1][0]
					if (bulgepos - blue.segment_dict[hairpin[1]].bp_data[0][0] - n_bulges) %2 == 0:
						n_bulges += 1
						bulge_to_add = bulgepos
			if n_bulges == 1:
				sthb = HbondsBulgedPairWithOffset(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue,bulge_position=bulge_to_add,str2_offset=1) ; cst_fileout.write(sthb)
				cst_fileout.write('\n')
			elif n_bulges == 0:
				sthb = HbondsRegularHairpin(strand1=hairpin[0],strand2=hairpin[1],blueprint=blue) ; cst_fileout.write(sthb)
				cst_fileout.write('\n')
			elif n_bulges > 1:
				print("Something wrong with bulge count!")

	corner = cornerConstraints(capType=cap,blueprint=blue); cst_fileout.write(corner)

	cst_fileout.close()

	blue.remodel_all()
	blue.bp_data[0] = [1, 'A', 'L', '.']
	blue.dump_blueprint(blueprint_name)
        # this where the glycine in the trp corner goes
	#os.remove(blueprint_name)

	#############################################
	#-----------------------
	# Write modified xml
	#-----------------------
#	xml_out = open('input.xml','w')
#	for line in xml_lines:
#		xml_out.write(line)

#	os.chdir('../')
