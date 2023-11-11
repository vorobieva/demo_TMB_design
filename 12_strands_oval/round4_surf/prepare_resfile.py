import pandas as pd
import numpy as np
import glob
import os
import re
from shutil import copyfile
from collections import Counter
import pyrosetta as py

py.init()

for pdb in glob.glob("*.pdb"):
	pose = py.pose_from_pdb(pdb)
	seq = pose.sequence()
	pdb_name = pdb.split('.')[0]
	resfile = pdb_name + "/surface.resfile"
	with open("surf.resfile",'r') as resfile_in:
		with open(resfile, 'w') as resfile_out:
			for line in resfile_in:
				if "#surface" in line:
					vals = line.split(" ")
					if seq[int(vals[0])-1] == "G":
						resfile_out.write("%s A PIKAA G\n" %(vals[0]))
					elif seq[int(vals[0])-1] == "S":
						resfile_out.write("%s A PIKAA S\n" %(vals[0]))
					elif seq[int(vals[0])-1] == "T":
						resfile_out.write("%s A PIKAA T\n" %(vals[0]))
					else:
						resfile_out.write(line)
				else:
					resfile_out.write(line)
