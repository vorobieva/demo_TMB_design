import glob
import os
import shutil

for pdb_file in glob.glob("*_input_*/*.pdb"):
	file_name = os.path.splitext(pdb_file)
	vals = file_name[0].split("_")
	num=vals[-1]
	num = num.lstrip("0")
	new_pdb = ("_").join(vals[0:-1])+"_"+num+".pdb"
	shutil.move(pdb_file,new_pdb)
