import glob
import itertools
import shutil
import re
import os

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

for pdb_file in glob.glob("*.pdb"):
	pdb_name = pdb_file.split(".")[0]
	nets = {}
	Y10 = []
	Y96 = []
	num_nets = len(glob.glob(pdb_name+"/*.cst"))
	for i in range (1,num_nets+1):
		net = pdb_name + "/" + pdb_name + "_0001_network_" + str(i)
		with open(net+".cst", 'r') as cst_file:
			net_name = ""
			cst = ""
			for line in cst_file:
				if "#network_" in line:
					vals = line.split()
					net_name = vals[2]
					if "Y_10" in net_name:
						Y10.append((net_name, i))
					elif "Y_96" in net_name:
						Y96.append((net_name, i))
				elif "#AtomPair" in line:
					cst = line[1:]
			nets[i] = (net_name, cst)
	prod = list(itertools.product(Y10, Y96))
	print(pdb_file, prod)
	for p in prod:
		n_comb = ",".join(n[0] for n in p)
		comb = [x[1] for x in p]
		c_comb = "".join(nets[y][1] for y in comb)
		nets[max(nets, key=int) + 1] = (n_comb, c_comb)
	
	for key in nets:
		if len(nets[key][0].split(",")) == 4: # Add this line if only combinations of motifs are to be considered (exclude single motifs)
			source_pdb = pdb_name + "/" + pdb_name + "_0001_network_" + str(key) + ".pdb"
			new_path = "../round2/"+pdb_name+"_"+str(key)+"/"
			try:
				os.mkdir(new_path)
			except OSError:
				print ("Creation of the directory %s failed" % new_path)	
	
			muts = nets[key][0].split(",")
			pairs = list(chunks(muts,2))
			m ={}
			for j in range(0, len(pairs)):
				m_name = ",".join(pairs[j])
				for k, d in nets.items():
					if nets[k][0] == m_name:
						s_pdb = pdb_name + "/" + pdb_name + "_0001_network_" + str(k) + ".pdb"
						res_id = re.findall("\d+", m_name)
						for res in res_id:
							coord = ""
							with open(s_pdb, 'r') as template:
								for line in template:
									if "ATOM" in line:
										if line.split()[5] == res:
											coord += line
							m[res] = coord
		
			with open(pdb_file, 'r') as template:
				with open(new_path+pdb_name+"_"+str(key) + ".pdb", 'w') as out_f:
					for l in template:
						if "ATOM" in l and l.split()[5] in m:
							if l.split()[2] == "N":
								out_f.write(m[l.split()[5]])
							else:
								pass
						elif "ATOM" in l:
							out_f.write(l)
					for r in m:
						out_f.write("REMARK PDBinfo-LABEL:   %s HBNet\n" %(r))

			with open(new_path+"cst", "w") as cst_file:
				for pair in pairs:
					p = sorted(pair)
					print(p)
					cst_file.write("Dihedral C %s CA %s CB %s CG %s CIRCULARHARMONIC -1.22 0.30\n" %(p[1].split("_")[2],p[1].split("_")[2],p[1].split("_")[2],p[1].split("_")[2]))
				cst_file.write(nets[key][1])
		
			motif_res = []
			for pair in pairs:
				motif_res.append(pair[0].split("_")[2])
				motif_res.append(pair[1].split("_")[2])

			with open("../all.resfile", "r") as template:
				with open(new_path+"core.resfile", "w") as resfile:
				# write a resfile with (a) beta-turn fixed sequences and (b) the sequences and positions of the found motifs.
					for l in template: 
						resfile.write(l)
						if l.startswith("start"):
							break
					for l in template:
						if l == "\n":
							resfile.write(l)
						else:
							if "#surface" in l:
								pass
							elif l.split()[0] in motif_res:
								pass
							else:
								resfile.write(l)
#					resfile.write("ALLAA\n\nstart\n\n2 A PIKAA NTEDPGRKQS\n6 A PIKAA STDN\n\n16 A PIKAA DENST\n17 A PIKAA AES\n18 A PIKAA D\n19 A PIKAA G\n\n32 A PIKAA STD\n33 A PIKAA P\n34 A PIKAA DEHTY\n\n44 A PIKAA DENST\n45 A PIKAA AES\n46 A PIKAA D\n47 A PIKAA G\n\n60 A PIKAA STD\n61 A PIKAA P\n62 A PIKAA Y\n\n72 A PIKAA DENST\n73 A PIKAA AES\n74 A PIKAA D\n75 A PIKAA G\n\n109 A POLAR\n\n90 A PIKAA STD\n91 A PIKAA P\n92 A PIKAA DEHTY\n\n104 A PIKAA DENST\n105 A PIKAA AES\n106 A PIKAA D\n107 A PIKAA G\n\n122 A PIKAA STD\n123 A PIKAA P\n124 A PIKAA DEHTY\n\n134 A PIKAA DENST\n135 A PIKAA AES\n136 A PIKAA D\n137 A PIKAA G\n\n150 A PIKAA STD\n151 A PIKAA P\n152 A PIKAA DEHTY\n\n162 A PIKAA DENST\n163 A PIKAA AES\n164 A PIKAA D\n165 A PIKAA G\n")
					for pair in pairs:
						resfile.write("%s A PIKAA %s\n%s A PIKAA %s\n" %(pair[0].split("_")[2],pair[0].split("_")[1],pair[1].split("_")[2],pair[1].split("_")[1]))	



