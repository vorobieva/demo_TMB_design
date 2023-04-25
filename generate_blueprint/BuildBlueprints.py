""
#Script for writing blueprint files from a string defining the protein topology
#H[n-m]L[n-m]
#Length: xxx
""
# USE: python build_blueprints.v2.py -xml template_bb+design.xml -blueresfile

import itertools
import re
import sys
import os
import copy
from Blueprint import Blueprint 
from itertools import chain
from argparse import ArgumentParser
import numpy as np

def Bulged(strand):
        flag=False ; bulgepos=None
        for i in range(1,len(strand.bp_data)-1):
                prev_res = strand.bp_data[i-1]
                res = strand.bp_data[i]
                next_res = strand.bp_data[i+1]
                if 'EA' == res[2] and 'EB' == prev_res[2] and 'EB' == next_res[2]:
                    flag=True
                    bulgepos = res[0]
        return bulgepos

def Kinked(strand):
        flag=False ; kinkpos=None
        for i in range(1,len(strand.bp_data)-1):
                prev_res = strand.bp_data[i-1]
                res = strand.bp_data[i]
                next_res = strand.bp_data[i+1]
                if 'EA' == res[2] and 'EB' == prev_res[2] and 'EB' == next_res[2]:
                    flag=True
                    bulgepos = res[0]
        return bulgepos

def MakePlainBlueprint(ss,comb,bluefile):
	c = comb     
	out_file = open(bluefile,'w')
	struct = {'H':'HA', 'E':'EB', 'L':'LA'}
	total_length=sum(c)
	k=0
	curr_ss = ss[k]
	for i in range(1,total_length+1):
		if i>sum(c[:k+1]):
			k+=1
			curr_ss = ss[k]
		out_file.write('0  V  %s  R\n' %(curr_ss))
	out_file.close()

#------------
def MakeRefBlueprint(ss,comb,**kwargs):
	refbluefile = kwargs.get('refblue')
	bulges = kwargs.get('bulges',None)
	strandABEGO = kwargs.get('strandsABEGO',None)
	mostCommonLoopABEGO =kwargs.get('mostCommonLoopABEGO', False)
	cap=kwargs.get('cap',None)
	ABEGOString = [''] * len(ss)
	if strandABEGO:
		for idx,segment in enumerate(ss):
			if segment == 'E':
				ABEGOString[idx] = 'B'
	c = comb     
	out_file = open(refbluefile,'w')
	struct = {'H':'HA', 'E':'EB', 'L':'LA'}
	total_length=sum(c)
	k=0
	curr_ss = ss[k]
	curr_abego = ABEGOString[k]
	for i in range(1,total_length+1):
		if i>sum(c[:k+1]):
			k+=1
			curr_ss = ss[k]
			curr_abego = ABEGOString[k]
		out_file.write('0  V  %s%s  R\n' %(curr_ss,curr_abego))
	out_file.close()

        # Put bulges
	blue = Blueprint(refbluefile)
	bluelist = []
	
	if bulges:
		for seg in blue.segments:
			if seg.id in bulges.keys():
				num_bulges = len(bulges[seg.id])
				for j, res in enumerate(seg.bp_data):
					for bulge in range(0, num_bulges):
						seg.bp_data[bulges[seg.id][bulge]-1][2] = 'EA'
#       	                	        	if j == bulges[seg.id][bulge]-1:
#       	                        	        	res[2] = 'EA'
	if mostCommonLoopABEGO:
		for seg_id, seg in enumerate(blue.segments):
			if seg.sstype == 'L':
			# segment is a beta-turn
				try:
					if blue.segments[seg_id-1].sstype == 'E' and blue.segments[seg_id+1].sstype == 'E':
						if len(seg.bp_data) == 3:
							turnABEGO = ['A','A','G']
							for line in range(len(seg.bp_data)):
								seg.bp_data[line][2] = "L" + turnABEGO[line]
						elif len(seg.bp_data) == 4:
							turnABEGO = ['A','A','A','G']
							for line in range(len(seg.bp_data)):
								seg.bp_data[line][2] = "L" + turnABEGO[line]
						elif len(seg.bp_data) == 5:
							turnABEGO = ['A','A','A','G','G']
							for line in range(len(seg.bp_data)):
								seg.bp_data[line][2] = "L" + turnABEGO[line]
						elif len(seg.bp_data) == 2:
							if blue.segments[seg_id-1].bp_data[-2][2] == 'EA':
								turnABEGO = ['A','A']
								for line in range(len(seg.bp_data)):
									seg.bp_data[line][2] = "L" + turnABEGO[line]
				except IndexError:
					continue 
	
	if cap: 
		if cap==2:
			blue.bp_data[3][2] = "LB"
			blue.bp_data[4][2] = "LG"

	blue.dump_blueprint(refbluefile)
	#os.chdir('../')

#---------------
def Shift(**kwargs):
	refbluefile = kwargs.get('refblue')
	seg1 = kwargs.get('seg1')
	seg2 = kwargs.get('seg2')
	refblue = Blueprint(refbluefile)
	nseg1 = len( refblue.segment_dict[seg1].bp_data )
	nseg2 = len( refblue.segment_dict[seg2].bp_data )
	shift = nseg1 - nseg2
	return shift
#---------------
def MakeFirstBlueprint(**kwargs):
	ctail = [[0, 'V', 'L', 'R']]
	refbluefile = kwargs.get('refblue')
	segments = kwargs.get('segments')
	newbluefile = kwargs.get('newblue')
	ss_pairing = kwargs.get('ss_pairing')
	hs_pairing = kwargs.get('hs_pairing')
	hh_pairing = kwargs.get('hh_pairing')
	seg_adapt = kwargs.get('adapt')
	ntail = kwargs.get('ntail',None)
	if not ntail:
		ntail = [[0, 'V', 'L', 'R']]

	refblue = Blueprint(refbluefile)

	if seg_adapt:
		seg1 = seg_adapt.keys()[0] # strand to adapt
		seg2 = seg_adapt[seg1] # referemce strand		
		nseg1 = len( refblue.segment_dict[seg1].bp_data )
		nseg2 = len( refblue.segment_dict[seg2].bp_data )
		shift = nseg1 - nseg2
	else:
		shift=0
	bp_data_new = []
	for seg in segments:
		if shift == 0:
			for res in refblue.segment_dict[seg].bp_data:
				bp_data_new.append(res)
		elif shift > 0:
			if seg==seg1:
				for k in range(0,nseg1-shift):
					bp_data_new.append(refblue.segment_dict[seg1].bp_data[k])
			else:
				for res in refblue.segment_dict[seg].bp_data:
					bp_data_new.append(res)
		elif shift < 0:
			if seg==seg2:
				for k in range(shift,nseg2):
					bp_data_new.append(refblue.segment_dict[seg2].bp_data[k])
					bp_data_new.append(refblue.segment_dict[seg1].bp_data[k])
			else:
				for res in refblue.segment_dict[seg].bp_data:
					bp_data_new.append(res)

	refblue.bp_data = ntail + bp_data_new + ctail

	# write strand pairing	
	ss_line = "SSPAIR "
	# Get first strand to set "start"
	segments.sort()
	npairs=0
	for i,key in enumerate(ss_pairing.keys()):
		for st in ss_pairing[key]:
			s1=key
			s2 = st[:st.find('.')]
			orient = st[-1]
			# compute effective length of strands for calculating register shift
			#------
			if Bulged(refblue.segment_dict[s1]):
				n1 = len( refblue.segment_dict[s1].bp_data ) - 1
			else:
				n1 = len( refblue.segment_dict[s1].bp_data )
			if Bulged(refblue.segment_dict[s2]):
				n2 = len( refblue.segment_dict[s2].bp_data ) - 1
			else:
				n2 = len( refblue.segment_dict[s2].bp_data )
			#------
			if orient == 'A':
				shift = n1-n2
		
			if npairs==0:
				ss_line += '%s-%s.%s.%s' %(int(s1[1:]),int(s2[1:]),orient,shift)	
			else:
				ss_line += ';%s-%s.%s.%s' %(int(s1[1:]),int(s2[1:]),orient,shift)	
				npairs+=1

	header = [ss_line]
	if hs_pairing != None:
		header.append( hs_pairing )

	if hh_pairing !=None:
		header.append( hh_pairing )

	refblue.dump_blueprint(newbluefile,header_lines=header)
	
	# Abego 'B' for building strand
	os.system("sed  's/ E / EB/' %s > %s.b" %(newbluefile,newbluefile))
				
			
#-------------------
def AddSegmentToBlueprint(**kwargs):
	tail = [[0, 'V', 'L', 'R']]
	refbluefile = kwargs.get('refblue')
	segments = kwargs.get('segments')
	blue0file= kwargs.get('blue0')
	newbluefile = kwargs.get('newblue')
	append = kwargs.get('append')
	insert_between_first = kwargs.get('insert_between_first')
	insert_between_last = kwargs.get('insert_between_last')
	ss_pairing = kwargs.get('ss_pairing')
	ss_pairing_shift = kwargs.get('ss_pairing_shift')
	hs_pairing = kwargs.get('hs_pairing')
	hh_pairing = kwargs.get('hh_pairing')
	seg_abego = kwargs.get('seg_abego')
	specific_abego = kwargs.get('specific_abego')
	insert = kwargs.get('insert')
	only_remodel = kwargs.get('only_remodel',False)

	blue0 = Blueprint(blue0file)
	blue0.reindex_blueprint(start=1)
	blue0.freeze_all()

	refblue = Blueprint(refbluefile)
	bp_data_new = []
	for seg in segments:
		for res in refblue.segment_dict[seg].bp_data:
			bp_data_new.append(res)		

	if append==True:
		blue0.bp_data[-2][3]='R'
		blue0.bp_data[-1][3]='R' ; blue0.bp_data[-1][2]='L'
		bp_data_new = blue0.bp_data + bp_data_new[1:] + tail

	elif append==False:
		blue0.bp_data[0][3]='R' #; blue0.bp_data[0][2]='L'
		blue0.bp_data[1][3]='R'
		if segments[-1][0] != 'L':
			blue0.bp_data[0][2]= bp_data_new[-1][2]
		bp_data_new =  tail + bp_data_new[:-1] + blue0.bp_data

	elif insert_between_first:
		# replace original residues by insert ones. so we can use original bp as blue0
		index1 = insert_between_first - 1
		index2 = insert_between_last - 1
		if index1>1:  # we need at least one point fixed (not R)
			blue0.bp_data[index1][3]='R' #; blue0.bp_data[index1][2]='LX'
		
		bp_data_new = blue0.bp_data[:index1+1] + bp_data_new # + blue0.bp_data[index2:] 
		shift = index2-index1-1
		for k in range(index2,len(blue0.bp_data)):
			blue0.bp_data[k][0]-=shift
		blue0.bp_data[index2][3]='R'
		bp_data_new = bp_data_new + blue0.bp_data[index2:]

	else:
		bp_data_new = blue0.bp_data		

	blue0_top = blue0.topology() # for abego conversion later
	blue0_cp = copy.deepcopy(blue0)

	# if we dont add residues then we dont need to replace bp_data. it is already update with R
	if only_remodel==False or isinstance(only_remodel,list):
		blue0.bp_data = bp_data_new

	# write strand pairing
	if ss_pairing != None:
		ss_line = "SSPAIR "
		keys = ss_pairing.keys()
		keys.sort()
		npairs=0
		for i,key in enumerate(keys):
			for st in ss_pairing[key]:
				s1=key
				s2 = st[:st.find('.')]
				orient = st[-1]

				# compute effective length of strands for calculating register shift
				#------
				if Bulged(refblue.segment_dict[s1]):
					n1 = len( refblue.segment_dict[s1].bp_data ) - 1
				else:
					n1 = len( refblue.segment_dict[s1].bp_data )
				if Bulged(refblue.segment_dict[s2]):
					n2 = len( refblue.segment_dict[s2].bp_data ) - 1
				else:
					n2 = len( refblue.segment_dict[s2].bp_data )
				#------

				#if orient == 'A':
				#shift = n1-n2
				shift=99	
				if npairs==0:
					ss_line += '%s-%s.%s.%s' %(int(s1[1:]),int(s2[1:]),orient,shift)
				else:
					ss_line += ';%s-%s.%s.%s' %(int(s1[1:]),int(s2[1:]),orient,shift)
				npairs+=1

	header=[]
	if ss_pairing_shift != None:
		header.append( ss_pairing_shift )
	else:
		header.append( ss_line )
	if hs_pairing != None:
		header.append( hs_pairing )

	if hh_pairing !=None:
		header.append( hh_pairing )

	blue0.dump_blueprint(newbluefile, header_lines=header)


	# insertion
	if insert:
		seg = insert.keys()[0]
		shift = insert[seg]
		blue = Blueprint('%s' %(newbluefile))
		blue_aux = copy.deepcopy(blue)
		#print blue.segment_dict[seg].bp_data
		blue_aux.reindex_blueprint()
		insert_index = blue_aux.segment_dict[seg].bp_data[-1][0] # we want to insert just before the 1st residue of the next seg, so that we insert at the very end

		for k in range(shift):
			blue.bp_data.insert(insert_index, [0, 'V', 'E', 'R'])
			#print blue.segment_dict[seg].bp_data

		# FOR NEW CONSTRAINTS WE NEED 7 RESIDUES TO USE BEND CONSTRAINTS (5 + 2 INSERTED)
		#for i in range(1,5+1):
		#blue.segment_dict[seg].bp_data[-i][3]='R'

		blue.dump_blueprint(newbluefile,header_lines=header)



	# abego loop
	# Adapt global abego motif to current stage
	#---------------------------------
	if seg_abego != None:
		new_abego={};conversor={}
		r=re.compile('[HEL]')
		top= blue0_cp.topology()
		curr_ss = r.findall(top)
		counter=0
		for seg in segments:
			ss = seg[0]
			if append:
				if ss=='L':
					newindex = curr_ss.count(ss)
				else:
					newindex = curr_ss.count(ss) + 1
					curr_ss.append(ss)
			elif append==False:
				curr_ss.insert(counter,ss)
				newindex = curr_ss[:counter+1].count(ss)
				if ss=='L': newindex+=1 # because Nter is loop L1 (and is not added through segments)
				counter+=1
			elif insert_between_first:	
				r=re.compile('[HEL]\d') ; seg_list = r.findall(top)
				flag=False
				for k,sg in enumerate(seg_list): # find segment whhere insert_first is located
					for res in blue0_cp.segment_dict[sg].bp_data:
						if res[0]==insert_between_first+1:
							flag=True
							break
					if flag:
						break
	
				first_seg=k
				newindex = curr_ss[:first_seg+1+counter].count(ss)
				counter+=1

			conversor[seg]='%s%s' %(ss,newindex)

		for seg in seg_abego:
			new_abego[conversor[seg]] = seg_abego[seg]

		
		blue = Blueprint('%s' %(newbluefile))
		if new_abego != None:
			for seg in new_abego.keys():
				abego=new_abego[seg]
				for i,res in enumerate(blue.segment_dict[seg].bp_data):
					res[2]+=abego[i]

		blue.dump_blueprint(newbluefile,header_lines=header)
	
	if specific_abego: # only abegos for specific positions
		blue = Blueprint('%s' %(newbluefile))
		for seg in specific_abego:
			pos,letter =  specific_abego[seg]
			blue.segment_dict[seg].bp_data[pos][2]+=letter

		blue.dump_blueprint(newbluefile,header_lines=header)

	# only remodel.
	blue = Blueprint('%s' %(newbluefile))
	if isinstance(only_remodel,list):
		for seg in only_remodel:
			blue.remodel_segment(id=seg)

		blue.dump_blueprint(newbluefile,header_lines=header)
		
	# Abego 'B' for building strand
	os.system("sed  's/ E / EB/' %s > %s.b" %(newbluefile,newbluefile))

#-------------------
def write_dummy_pdb(filename):
    dummy_pdb = 'ATOM      1  N   GLY A   1       0.346   1.210   0.744  1.00  0.00 \n\
ATOM      2  CA  GLY A   1       1.687   1.135   0.174  1.00  0.00 \n\
ATOM      3  C   GLY A   1       2.383   2.488   0.222  1.00  0.00 \n\
ATOM      4  O   GLY A   1       2.996   2.918  -0.752  1.00  0.00 \n\
ATOM      5 1H   GLY A   1      -0.448   0.981   0.185  1.00  0.00 \n\
ATOM      6 2H   GLY A   1       0.109   0.647   1.536  1.00  0.00 \n\
ATOM      7 3H   GLY A   1       0.005   2.079   1.101  1.00  0.00 \n\
ATOM      8 1HA  GLY A   1       2.277   0.413   0.741  1.00  0.00 \n\
ATOM      9 2HA  GLY A   1       1.615   0.810  -0.864  1.00  0.00 \n\
ATOM     10  N   GLY A  2       2.284   3.158   1.368  1.00  0.00 \n\
ATOM     12  CA  GLY A  2       3.918   4.450   2.676  1.00  0.00 \n\
ATOM     13  C   GLY A  2       3.551   4.379   3.850  1.00  0.00 \n\
ATOM     14  O   GLY A  2       1.859   5.564   1.752  1.00  0.00 \n\
ATOM     15 1H   GLY A  2       1.061   5.930   0.512  1.00  0.00 \n\
ATOM     16 2H   GLY A  2      -0.012   6.930   0.747  1.00  0.00 \n\
ATOM     17 3H   GLY A  2      -0.863   7.182  -0.405  1.00  0.00 \n\
ATOM     18 1HA  GLY A  2      -0.564   8.037  -1.404  1.00  0.00 \n\
ATOM     19 2HA  GLY A  2       0.540   8.749  -1.379  1.00  0.00 \n'
    out = open(filename, 'w')
    out.write(dummy_pdb)

def Bulges(blue):
        bulge_list=[]
        for i,res in enumerate(blue.bp_data):
                if res[2]=='EA':
                        bulge_list.append(res[0])
        return bulge_list

def XMLReplaceTagsValues(**kwargs):
        xml_lines = kwargs.get('xml_lines')
        identifier = kwargs.get('identifier')
        tags = kwargs.get('tags')
        values = kwargs.get('values')

        for i,line in enumerate(xml_lines):
                if identifier in line:
                        for tag,value in zip(tags,values):
                                line = line.replace('%s' %(tag),'%s' %(value))
                                xml_lines[i] = line

def XMLReplaceXXXYYY(**kwargs):
	xml_lines = kwargs.get('xml_lines')
	identifier = kwargs.get('identifier')
	xxx = kwargs.get('xxx')
	yyy = kwargs.get('yyy')

	for i,line in enumerate(xml_lines):
		if identifier in line:
			if xxx:
				line = line.replace('xxx','%s' %(xxx))
				xml_lines[i] = line

			if yyy:
				line = line.replace('yyy','%s' %(yyy))
				xml_lines[i] = line


def AmbiguousConstraints(list1,list2):
	st='AmbiguousConstraint\n'
	for res1 in list1:
		for res2 in list2:
			st += "AtomPair N %i O %i BOUNDED 3.5 4.5 0.5\n" %(res1,res2)
	st+='END_AMBIGUOUS\n'	
	return st

def ReplaceLine(line,word1,word2):
	if word1 in line:
		line.replace(word1,word2)


# super helpful utility functions to handle the parsing of the beta-barrel topology string and different sselements lengths. code taken from 
# https://stackoverflow.com/questions/5704931/parse-string-of-integer-sets-with-intervals-to-list

def group_to_range(group, strand=False):
    group = ''.join(group.split())
    sign, g = ('-', group[1:]) if group.startswith('-') else ('', group)
    r = g.split('-', 1)
    r[0] = sign + r[0]
    r = sorted(int(__) for __ in r)
    increment = 2 if strand else 1
    return range(r[0], 1 + r[-1], increment)

def rangeexpand(txt, strand=False):
    txt = txt.replace('[','')
    txt = txt.replace(']','')
    ranges = chain.from_iterable(group_to_range(__, strand) for __ in txt.split(','))
    return sorted(set(ranges))

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def getSpecificLengths(topol):
    """
    This function convert the mix of discrete lengths and ranges specified for a given secondary structure element in the topol string into discrete lengths only.
    Args:
	topol (string): a user-specified string describing the sequence of ssecondary structure elements and lengths to consider. 
    Returns:
	ss (list): a list of the secondary structure elements described as a 1-letter code ['L','E',...]
	itemizedLength (list): Each element of the list is a list of integers corresponding to all possible lengths for a given secondary structure element. 
    """
    elements = re.compile("[HEL]")
    ss = elements.findall(topol)
    lengths = elements.split(topol)[1:]
    itemizedLength = [rangeexpand(length, ss_elem == "E") for ss_elem, length in zip(ss, lengths)] 
    return ss, itemizedLength

def identify_elongation_targets(aa, chunk_size=1):
    """
    For a given type of secondary structure, this function finds the elements that have been specified different lenghts and that will be elongated. 
    It treats the strands as pairs that has to be elongated together.
    Args:
        aa (list): a list of tuples of the form ('L', [3, 4, 5]) to specify the ssElement
            type and legal lengths for it.
    
    Returns:
	elongated_strand_pairs (list): Each element of the list is a tuple of integer(s) that corresponds to the index of the secondary structure element(s) to elongate. 
    """
    elongated_strand_pairs = []
    for i, chunked_ss in enumerate(chunks(aa, chunk_size)):
        assert(len({len(ss[-1]) for ss in chunked_ss}) == 1)
        if len(chunked_ss[0][-1]) > 1:
            elongated_strand_pairs.append(tuple(i * chunk_size + j for j in range(chunk_size)))
    return elongated_strand_pairs

def elongate_sec_struc(ssElemLeng, elongation_targets):
    """
    Examples:
    
    Args:
        ssElemLeng (list): a list of tuples of the form ('E', [len1, len2,..., lenN])
        elongation_targets (list): a list of indices (or groups of indices) to iterate over
    """
    elements = set()
    for elongation_target in elongation_targets:
        number_of_lengths = len(ssElemLeng[elongation_target[0]][-1])
        for i in range(number_of_lengths):
            elemLengths = []
            for sNo, (ssElem, ssLength) in enumerate(ssElemLeng):
                if sNo in elongation_target:
                    elemLengths.append(ssLength[i])
                    continue
                elemLengths.append(ssLength[0])
            elements.add(tuple(elemLengths))
    return list(elements)

def getShortestIdealBarrelFromInput(ss, itemizedLength):
    """This function assumes that the starting template is with the shorter hairpins and loops.
    Elongations are additional to the idealized barrel
    """
    return [(ssElem, min(itemLength)) for ssElem, itemLength in zip(ss, itemizedLength)]

def calculate_kink_combinations(input_kink_positions,ordered_cbeta_strips, max_kinks, min_kinks):
	possible_kink_intervals = input_kink_positions

	all_to_all = list(itertools.product(possible_kink_intervals, repeat=len(ordered_cbeta_strips)))
	filtered_all_to_all = []
	for kink_combination in all_to_all:
    		if (len(list(chain.from_iterable(kink_combination))) <= max_kinks) and (len(list(chain.from_iterable(kink_combination))) >= min_kinks):
        		filtered_all_to_all.append(kink_combination)

	return filtered_all_to_all

def findDoubleOffsets(nStrand, shearNo, ideal_barrel):
    	"""this function finds the double offsets between the strands based on the ideal configuraion
    	of the barrel.
    
    	Notes:
           This function assumes the target barrel is ideal 
           with ideal strand offsets. Strands that require a beta bulge will have a residue added
           to them later to account for this.
    	"""
    
    	strands_per_hairpin = 2
    	strands = [leng for ss_elem, leng in ideal_barrel if ss_elem == "E"]
    	shear_no = sum([b - a for a, b in chunks(strands, strands_per_hairpin)])

    	assert(len(strands) == nStrand)
    	assert(shear_no == shearNo)
    
    	double_offset_strand_idx = [(i * strands_per_hairpin) + 1 for i, (a, b) in 
                                 enumerate(chunks(strands, strands_per_hairpin)) if b - a == 4]
    	return double_offset_strand_idx

def findCbetaStrips(nStrand, shearNo, topol):
	# The first few lines are there to verify that the specified topology matches the target type of barrel specified by the user.
	ss, itemizedlength = getSpecificLengths(topol)
	ideal_barrel = getShortestIdealBarrelFromInput(ss, itemizedlength)
	dosStrands = findDoubleOffsets(nStrand, shearNo, ideal_barrel)
	assert(len(dosStrands) == (shearNo - nStrand) / 2)

	_aa = [strandDir * shear for shear in range(2, shearNo + 1, 2) for strandDir in (-1, 1)]
	cbeta_strip_offsets_unfiltered = [[_aa[strand - shear] for strand in range(len(_aa))] 
                                      for shear in range(0, shearNo, 2)]

	cbetaStrips = []
	skipPair = [int(i/2) for i in sorted(dosStrands)]
	for strip in cbeta_strip_offsets_unfiltered:
		cbeta_strip_offsets_filtered = [offsets for hairpin, offsets in enumerate(chunks(strip, 2)) if hairpin not in skipPair]
		cbetaStrips.append(list(chain.from_iterable(cbeta_strip_offsets_filtered)))
	
	cbetaStrips = [bb[-1:] + bb[:-1] for bb in cbetaStrips]
    	
	# reorder Cbeta-strips and add strand numbers to make it easier to fetch positions
	orderedCbetaStrip = []

	for strip in cbetaStrips:
		start_strand = strip.index(max(strip))+1
		pos_in_strand = []
		for pos in enumerate(strip):
			pos_in_strand.append(pos)
		pos_in_strip = pos_in_strand[start_strand:]+pos_in_strand[:start_strand]
		orderedCbetaStrip.append(pos_in_strip)
	return orderedCbetaStrip

def placeBulgesOnElongatedStrands(complete_strand_lengths, bulges, strand_elongation_pairs):
	complete_strand_lengths_with_bulges = complete_strand_lengths
	for strand_1, strand_2 in strand_elongation_pairs:
    		bulge = {'E2':[-2],'E4':[-2],'E6':[-2]}
    		target_combo = max(complete_strand_lengths, key=lambda item: item[strand_2])
    		shortest_combo = min(complete_strand_lengths, key=lambda item: item[strand_2])
    		bulge_pos = target_combo[strand_2] - shortest_combo[strand_2] +1
    		target_combo_with_bulge = target_combo[:strand_2] + (target_combo[strand_2]+1,) + target_combo[strand_2+1:]
    		strand_idx = "E" + str(strand_2+1)
    		if strand_idx in bulge:
        		bulge[strand_idx] = bulge[strand_idx] + [bulge_pos]
    		else:
        		bulge[strand_idx] = [bulge_pos]
    		complete_strand_lengths_with_bulges.append(target_combo_with_bulge)
    		bulges.append(bulge)

	for comb, bulges_comb in enumerate(bulges):
    		for key in bulges_comb.keys():
        		for b_idx, b in enumerate(bulges_comb[key]):
            			if b < 0 :
                			strand_idx = int(key[1])-1
                			bulges_comb[key][b_idx] = complete_strand_lengths_with_bulges[comb][strand_idx] + b + 1

	return complete_strand_lengths_with_bulges, bulges

def recalculateBulgePos(bulges, complete_strand_lengths):
	for comb, bulges_comb in enumerate(bulges):
		for key in bulges_comb.keys():
			for b_idx, b in enumerate(bulges_comb[key]):
				if b < 0 :
					strand_idx = int(key[1])-1
					bulges_comb[key][b_idx] = complete_strand_lengths[comb][strand_idx] + b + 1
	return bulges
    		
def add_tryptophan_corner(ss, ss_leng_comb):
	ssWCorner = ['L','H','L']
	ss_with_WCorner = ssWCorner + ss
	lengWCorner = (2, 4, 2)
	# make a tuple of lengths of each ss element -- subtract one residue from the
	# first and last strands (indices 0 and -2 in the length_combo tuple) to account
	# for the trp corner
	leng_with_WCorner = tuple(lengWCorner + (ss_leng_comb[0] - 1,) + ss_leng_comb[1:-2] + (ss_leng_comb[-2] - 2,) + (ss_leng_comb[-1],))
	return ss_with_WCorner, leng_with_WCorner

def add_bserine_corner(ss, ss_leng_comb):
        ssSCorner = ['L']
        ss_with_SCorner = ssSCorner + ss
        lengSCorner = (5,)
        # make a tuple of lengths of each ss element -- add one residue (the serine) to the
        # first strand 
	# (ss index 0 in the length_combo tuple) to account
        leng_with_SCorner = tuple(lengSCorner + ss_leng_comb[:])
        return ss_with_SCorner, leng_with_SCorner

def getCombinationsForBetaBarrels(topol,common_bulges):
	"""
	This finction calculates all possible combinations of lengths for beta-strands and loops, in a fashion specific to beta-barrels. 
	This means, beta-strands are elongated by pairs, since we assume antiparallel strands. 
	The beta-bulges positions are not specifed by the user but hard-coded and derived from the topology.
	Args:
		topol (str): a string describing the topology in terms of the secondary structure slements and possible lengths, specified by the user.
	Results:
	
	"""
	from copy import deepcopy

	ss, leng = getSpecificLengths(topol)
	bulges = common_bulges # {'E2':[-2],'E4':[-2],'E6':[-2]} #bulges positions on bottom hairpins

	ee = [(ss_elem, ssLength) for ss_elem, ssLength in zip(ss, leng) if ss_elem == "E"]
	ll = [(ss_elem, ssLength) for ss_elem, ssLength in zip(ss, leng) if ss_elem == "L"]

	strand_elongation_pairs = identify_elongation_targets(ee, chunk_size=2)
	loop_elongations = identify_elongation_targets(ll)

	if (len(e[1]) == 1 for e in ee) and (len(l[1]) == 1 for l in ll):
		# special case, no strand elngation required
		complete_strand_lengths = [tuple(e[1][0] for e in ee)]
		complete_loop_lengths = [tuple(l[1][0] for l in ll)]
	elif (len(e[1]) == 1 for e in ee):
		complete_strand_lengths = [tuple(e[1][0] for e in ee)]
		complete_loop_lengths = elongate_sec_struc(ll, loop_elongations)
	elif (len(l[1]) == 1 for l in ll):
		complete_strand_lengths = elongate_sec_struc(ee, strand_elongation_pairs)
		complete_loop_lengths = [tuple(l[1][0] for l in ll)]
	else:
		complete_strand_lengths = elongate_sec_struc(ee, strand_elongation_pairs)
		complete_loop_lengths = elongate_sec_struc(ll, loop_elongations)

	# Now we add the bulges to all bottom hairpins (strands with index 1,3 and 5) and to elongated strands.
	bulges_comb = []
	for i in range(len(complete_strand_lengths)):
		bulges_comb.append(deepcopy(bulges))
		for key in bulges:
			strNo = int(key[1:]) - 1
			nextStr = int(key[1:])
			complete_strand_lengths[i] = complete_strand_lengths[i][:strNo] + (complete_strand_lengths[i][strNo]+1,) + complete_strand_lengths[i][nextStr:]
#			complete_strand_lengths[i] = complete_strand_lengths[i][:1] + (complete_strand_lengths[i][1]+1,) + complete_strand_lengths[i][2:]
#                	complete_strand_lengths[i] = complete_strand_lengths[i][:3] + (complete_strand_lengths[i][3]+1,) + complete_strand_lengths[i][4:]
#               	complete_strand_lengths[i] = complete_strand_lengths[i][:5] + (complete_strand_lengths[i][5]+1,) + complete_strand_lengths[i][6:]

	if len(strand_elongation_pairs) > 0 :
	# if some strand pairs have to be elongated, then we might have to add bulges
		complete_strand_lengths_with_bulges, complete_bulges_comb = placeBulgesOnElongatedStrands(complete_strand_lengths, bulges_comb, strand_elongation_pairs)
	elif len(strand_elongation_pairs) == 0:
		complete_strand_lengths_with_bulges, complete_bulges_comb = complete_strand_lengths, recalculateBulgePos(bulges_comb,complete_strand_lengths)

	combinations = []
	for strands, loops in itertools.product(complete_strand_lengths_with_bulges, complete_loop_lengths):
		combinations.append(tuple(chain.from_iterable(zip(strands, loops))))
	        
	resolved_bulges = [bulge for bulge in complete_bulges_comb for _ in range(len(complete_loop_lengths))]
	assert len(resolved_bulges) == len(combinations)
	return ss, combinations, resolved_bulges

def GetCombinations(topol):
	elements = re.compile("[HEL]")
	ss = elements.findall(topol)
	relengths = re.compile("(\d+)-(\d+)")
	relengths2 = re.compile("(\d+),(\d+)") # for specific lengths, not ranges
	lengths = relengths.findall(topol)
	lengths2 = relengths2.findall(topol)

	## for interpreting specific lengths
	index=-1
	specific=[]
	for st in topol:
		if st=='[':
			index+=1
		if st==',':
			specific.append(index)

	comb=[] ; j=0 ; k=0
	for i in range(len(ss)):

		if i in specific:
			frag = lengths2[j]
			comb.append( [int(l) for l in frag] )		
			j+=1
		else:
			fragment = lengths[k]
			comb.append(range(int(fragment[0]),int(fragment[1])+1))
			k+=1

	combinations = list(itertools.product(*comb))
	print('Number of combinations: %s' %(len(combinations)))

	return ss,combinations

def cornerConstraints(**kwargs):
	capType = kwargs.get("capType", 1)
	blue = kwargs.get("blueprint")

	if capType == 1:
		Turn2RefPos = blue.segment_dict["L4"].bp_data[0][0]
		Turn4RefPos = blue.segment_dict["L6"].bp_data[0][0]
		Turn6RefPos = blue.segment_dict["L8"].bp_data[0][0]
		ArgCornerRefPos = blue.segment_dict["E8"].bp_data[-2][0]
		CornerCst = "Dihedral N 8 CA 8 C 8 N 9 CIRCULARHARMONIC 2.35 0.25\nDihedral C 7 N 8 CA 8 C 8 CIRCULARHARMONIC 5.20 0.25\nDihedral N 7 CA 7 C 7 N 8 CIRCULARHARMONIC 5.75 0.25\nDihedral C 6 N 7 CA 7 C 7 CIRCULARHARMONIC 4.90 0.25\nAtomPair CA 11 O 8 BOUNDED 6.9 7.5 0.5\nAtomPair CA 7 CA %s BOUNDED 9.5 10.5 0.5\nAtomPair CA 7 CA %s BOUNDED 7.0 8.0 0.5\nAtomPair CA 6 CA %s HARMONIC 12.0 0.5\nAtomPair CA 6 CA %s HARMONIC 9.5 0.5\nAtomPair CA 6 CA %s HARMONIC 9.5 0.5\nAtomPair CA 7 CA %s BOUNDED 8.5 10.0 0.5\nAtomPair O 6 CA %s HARMONIC 8.5 0.5\n" %(Turn4RefPos,Turn6RefPos,Turn2RefPos,Turn4RefPos,Turn6RefPos,ArgCornerRefPos,ArgCornerRefPos)
	elif capType == 2:
		neighborBulgePos = blue.segment_dict["E2"].bp_data[-3][0]
		CornerCst = "AtomPair N 5 O %s HARMONIC 3.0 0.5\nAngle N 5 H 5 O %s CIRCULARHARMONIC 3.1 0.3\nAtomPair N 6 O 3 HARMONIC 3.0 0.5\nAngle N 6 H 6 O 3 CIRCULARHARMONIC 3.1 0.3\nDihedral C 2 N 3 CA 3 C 3 CIRCULARHARMONIC 3.65 0.25\nDihedral N 2 CA 2 C 2 N 3 CIRCULARHARMONIC 2.6 0.25\nDihedral C 1 N 2 CA 2 C 2 CIRCULARHARMONIC 5.25 0.25\n" %(neighborBulgePos,neighborBulgePos)
	else:
		CornerCst = ""
	return CornerCst

def HBondConstraints(donor,acceptor): # return string for cst file
    hb_ang = np.deg2rad(160.)
    hb_ang_tol=np.deg2rad(20.0)
    st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(donor,acceptor) 
    st+= "Angle N %i H %i O %i HARMONIC %3.1f %3.1f\n" %(donor,donor,acceptor,hb_ang,hb_ang_tol)
 #  st+= "Angle H %i O %i C %i HARMONIC %3.1f %3.1f\n" %(donor,acceptor,acceptor,hb_ang,hb_ang_tol)
    return st	

def HBondConstraintsBarrel(donor,acceptor): # return string for cst file
    hb_ang = np.deg2rad(160.)
    hb_ang_tol=np.deg2rad(20.0)
    st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(donor,acceptor)
    st+= "Angle N %i H %i O %i HARMONIC %3.1f %3.1f\n" %(donor,donor,acceptor,hb_ang,hb_ang_tol)
    return st

def CircularHBondConstraints(donor,acceptor): # return string for cst file
    hb_ang = np.deg2rad(180.)
    hb_ang_tol=np.deg2rad(20.0)
    st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(donor,acceptor)
    st+= "Angle N %i H %i O %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,donor,acceptor,hb_ang,hb_ang_tol)
#    st+= "Angle H %i O %i C %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,acceptor,acceptor,hb_ang,hb_ang_tol)
    return st

def CircularHBondConstraintsBarrel(donor,acceptor): # return string for cst file
    hb_ang = np.deg2rad(180.)
    hb_ang_tol=np.deg2rad(20.0)
    st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(donor,acceptor)
    st+= "Angle N %i H %i O %i CIRCULARHARMONIC %3.1f %3.1f\n" %(donor,donor,acceptor,hb_ang,hb_ang_tol)
    return st

def PairConstraints(a,b,value,tol,tag): # return string for cst file
    st = "AtomPair CA %i CA %i BOUNDED %3.1f %3.1f %3.1f 0.5 %s\n" %(a,b,value-tol,value+tol,tol/2,tag) 
    return st

def HarmonicPairConstraints(a,b,value,sd): # return string for cst file
    st = "AtomPair CA %i CA %i HARMONIC %3.1f %3.1f\n" %(a,b,value,sd)
    return st


def AngleConstraints(a,b,c,value,tol,tag): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Angle CA %i CA %i CA %i BOUNDED %3.1f %3.1f 0.5 %s\n" %(a,b,c,ang-ang_tol,ang+ang_tol,tag)
    return st

def HarmonicAngleConstraints(a,b,c,value,tol): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Angle CA %i CA %i CA %i HARMONIC %3.1f %3.1f\n" %(a,b,c,ang,ang_tol)
    return st


def CstTypeAngleConstraints(a,b,c,value,tol,cst_type): # return string for cst file
	st=''	
	if cst_type=='harmonic':
		ang = np.deg2rad(value)
		ang_tol = np.deg2rad(tol)
		st = "Angle CA %i CA %i CA %i HARMONIC %3.1f %3.1f\n" %(a,b,c,ang,ang_tol)
	elif cst_type=='bounded':
		ang = np.deg2rad(value)
		ang_tol = np.deg2rad(tol)
		sd = ang_tol/2
		tag="ang_%i.%i.%i" %(a,b,c)
		st = "Angle CA %i CA %i CA %i BOUNDED %3.2f %3.2f %3.2f 0.5 %s\n" %(a,b,c,ang-ang_tol,ang+ang_tol,sd,tag)
	return st


def DihedralConstraints(a,b,value,tol,tag): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)    
    sd = ang_tol/2
    st = "Dihedral CB %i CA %i CA %i CB %i BOUNDED %3.2f %3.2f %3.2f 0.5 %s\n" %(a,a,b,b,ang-ang_tol,ang+ang_tol,sd,tag)
    return st

def CaDihedralConstraints(a,b,c,d,value,tol): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    sd = ang_tol/2
    tag="ang_%i.%i.%i.%i" %(a,b,c,d)
    st = "Dihedral CA %i CA %i CA %i CA %i BOUNDED %3.2f %3.2f %3.2f 0.5 %s\n" %(a,b,c,d,ang-ang_tol,ang+ang_tol,sd,tag)
    return st

def CircularHarmonicCaDihedralConstraints(a,b,c,d,value,tol): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Dihedral CA %i CA %i CA %i CA %i CIRCULARHARMONIC %3.1f %3.1f\n" %(a,b,c,d,ang,ang_tol)
    return st

def HarmonicDihedralConstraints(a,b,value,tol): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Dihedral CB %i CA %i CA %i CB %i HARMONIC %3.1f %3.1f\n" %(a,a,b,b,ang,ang_tol)
    return st

def CircularHarmonicDihedralConstraints(a,b,value,tol): # return string for cst file
    ang = np.deg2rad(value)
    ang_tol = np.deg2rad(tol)
    st = "Dihedral CB %i CA %i CA %i CB %i CIRCULARHARMONIC %3.1f %3.1f\n" %(a,a,b,b,ang,ang_tol)
    return st


def AmbiguousCst(cst_lst):
        header = 'AmbiguousConstraint\n'
        for cst in cst_lst:
                header += cst
        header += 'END_AMBIGUOUS\n'
        return header


def MultiCst(cst_lst):
        header = 'MultiConstraint\n'
        for cst in cst_lst:
                header += cst
        header += 'END_MULTI\n'
        return header

def ConstraintsStrandCurvature(**kwargs):
	segment = kwargs.get("strand")
	positions = kwargs.get("positions")
	bend = float( kwargs.get("bend") )
	bend_tol = float( kwargs.get("bend_tol") )
	bend_bulge = float( kwargs.get("bend_bulge") )
	twist = float( kwargs.get("twist") )
	twist_tol = float( kwargs.get("twist_tol") )
	bluefile = kwargs.get("bluefile")

	blue = Blueprint(bluefile)
	blue.reindex_blueprint(start=1)
	seg = blue.segment_dict[segment]
	blue = Blueprint(bluefile)
	cst_st=''
	# bending
	if bend:
		if positions == None:
			positions = range( 2,len(seg.bp_data)-2)
		for i in positions:
			pos = seg.bp_data[i][0]
			if seg.bp_data[i][2] == 'EA': # bulge
				st = AngleConstraints(pos-2,pos,pos+2,180-bend_bulge,bend_tol,"bend%s.%i" %(segment,pos))
			else: # non-bulged
				st = AngleConstraints(pos-2,pos,pos+2,180-bend,bend_tol,"bend%s.%i" %(segment,pos))
			cst_st += st
	  
		pos1 = seg.bp_data[0][0]
		pos2 = seg.bp_data[-1][0]
		if len(seg.bp_data) % 2 ==0:
			cen1 = pos1 + len(seg.bp_data)/2
			cen2 = pos1 + len(seg.bp_data)/2 + 1
			st = AngleConstraints(pos1,cen1,pos2,180-bend_bulge,bend_tol,"edgebend%s.%i" %(segment,cen1))
			st = AngleConstraints(pos1,cen2,pos2,180-bend_bulge,bend_tol,"edgebend%s.%i" %(segment,cen2))
		else:
			cen = pos1 + len(seg.bp_data)/2
			st = AngleConstraints(pos1,cen,pos2,180-bend_bulge,bend_tol,"edgebend%s.%i" %(segment,cen))
		cst_st += st

		st = AngleConstraints(pos-2,pos,pos+2,180-bend_bulge,bend_tol,"bend%s.%i" %(segment,pos))
	# twisting
	if twist:
		if positions == None:
			positions = range( len(seg.bp_data)-2)
		for i in positions:
			pos1 = seg.bp_data[i][0]
			pos2 = pos1+2
			st = DihedralConstraints(pos1,pos2,twist,twist_tol,'dih%s.%i' %(segment,pos1))
			cst_st += st
	return cst_st


def RegularStrandCurvature(**kwargs):
	segment = kwargs.get("strand")

	level = kwargs.get("level") # 1 or 2

	bend = kwargs.get("global_bend",None)
	bend_tol = kwargs.get("global_bend_tol",10.0) 
	twist = kwargs.get("global_twist",None )
	twist_tol = kwargs.get("global_twist_tol",5.0 )

	bend_area = kwargs.get("bend_area_value",None )
	bend_area_value = kwargs.get("bend_area_value",None )
	bend_area_value_tol = kwargs.get("bend_area_value_tol",5.0)  # n-ter, middle, c-ter

	twist_area = kwargs.get("twist_area")  # n-ter, middle, c-ter
	twist_area_value = kwargs.get("twist_area_value",None )
	twist_area_value_tol = kwargs.get("twist_area_value_tol",5.0 )

	bend_positions = kwargs.get("bend_positions",None )
	bend_positions_value = kwargs.get("bend_positions_value",None )
	bend_positions_value_tol = kwargs.get("bend_positions_value_tol",None )

	constraint_type = kwargs.get("constraint_type","harmonic" )

	blue = kwargs.get("blueprint")
	blue.reindex_blueprint(start=1)
	seg = blue.segment_dict[segment]
	cst_st=''
	#################
	n = len(seg.bp_data)        

	############
	# BENDING
	############

	#----------------------
	# Set all triads for bend calculation
	#----------------------
	bend_triads = []
	step=level*2
	for i in range(0,n-step):
		if i+step*2 < n:
			pos1 = i
			pos2 = i + step*1
			pos3 = i + step*2
			bend_triads.append([pos1,pos2,pos3])

	#----------------------
	# Define bend areas
	#---------------------
	# By positions. This is used especially for positions paired to bulges, where we expect higher bend than the rest of triads
	bend_positions_triads=[]
	if bend_positions:
		for relpos in bend_positions:
			if relpos < 0: # when giving rel position from C-terminal (for E3)
				pos = n+relpos
			else:
				pos=relpos
			for triad in bend_triads:
				if pos==triad[1]: # central position of triad
					bend_positions_triads.append( triad )

	# By areas					
	bend_area_triads=[]
	if bend_area =='n-term':
		bend_area_triads = [ bend_triads[0] ]
	elif bend_area =='c-term':
		bend_area_triads = [ bend_triads[-1] ]
	elif bend_area =='center':
		index =  len(bend_triads)/2 - 1
		if len(bend_triads) % 2 == 0:
			bend_area_triads = [ pair for k in bend_triads[index:index+2] ]
		else:
			bend_area_triads = [ bend_triads[index] ]

	############
	# TWIST
	############

	#----------------------
	# Set all pairs for twist calculations
	#----------------------
	twist_pairs = []
	step=level*2
	for i in range(0,n-step):
		if i+step < n:
			pos1 = i
			pos2 = i + step*1
			twist_pairs.append([pos1,pos2])       

	#----------------------
	# Define twist areas
	#----------------------
	twist_area_pairs=[]
	if twist_area =='n-term':
		twist_area_pairs = [ twist_pairs[0] ]
	elif twist_area =='c-term':
		twist_area_pairs = [ twist_pairs[-1] ]
	elif twist_area =='center':
		index =  len(twist_pairs)/2 - 1
		if len(twist_pairs) % 2 == 0:
			twist_area_pairs = [ pair for pair in twist_pairs[index:index+2] ]
		else:
			twist_area_pairs = [ twist_pairs[index] ]
				
	#########################
	# INTRODUCE CONSTRAINTS
	#########################
	# After indentifying combination of positions for bending and twist... Put constraints
	# Calculate Bends
	for triad in bend_triads:
		a,b,c = triad	
		pos1 = seg.bp_data[b][0]
		pos2 = seg.bp_data[a][0]
		pos3 = seg.bp_data[c][0]         
		if bend_positions and triad in bend_positions_triads:
			st = CstTypeAngleConstraints(pos2,pos1,pos3,180-bend_positions_value,bend_positions_value_tol,constraint_type) ; cst_st += st
		elif bend_area_value and triad in bend_area_value_triads and triad not in bend_positions_triads:
			st = CstTypeAngleConstraints(pos2,pos1,pos3,180-bend_area_value,bend_area_value_tol,constraint_type) ; cst_st += st
		elif bend:
			st = CstTypeAngleConstraints(pos2,pos1,pos3,180-bend,bend_tol,constraint_type) ; cst_st += st

	# Calculate Twists
	for pair in twist_pairs:
		a,b = pair
		pos1 = seg.bp_data[a][0]
		pos2 = seg.bp_data[b][0]
		if twist_area_value and pair in twist_area_pairs:
			st = DihedralConstraints(pos1,pos2,twist_area_value,twist_area_value_tol,'dih%i.%i' %(pos1,pos2)) ; cst_st += st
		elif twist:
			st = DihedralConstraints(pos1,pos2,twist,twist_tol,'dih%i.%i' %(pos1,pos2)) ; cst_st += st 

	return cst_st

def BulgedStrandCurvature(**kwargs):
	segment = kwargs.get("strand")
	bend1 = kwargs.get("bend1",None)
	bend1_tol = kwargs.get("bend1_tol",5.0) 
	bend2 = kwargs.get("bend2",None)
	bend2_tol = kwargs.get("bend2_tol",5.0)
	blue = kwargs.get("blueprint")
	constraint_type = kwargs.get("constraint_type","harmonic" )

	seg = blue.segment_dict[segment]
	nres = len(seg.bp_data)
	cst_st=''

	pos1 = Bulged(seg)
	pos2 = pos1+1
	bulgepos_fromN = pos1-seg.bp_data[0][0]+1
	bulgepos_fromC = pos1-seg.bp_data[-1][0]-1
	

	# level 1
	if bend1:
		if nres > 6 and ( bulgepos_fromN >=3 and bulgepos_fromC <=-4 ):
			st = CstTypeAngleConstraints(pos1-2,pos1,pos1+3,180-bend1,bend1_tol,constraint_type) ; cst_st += st
			st = CstTypeAngleConstraints(pos2-3,pos2,pos2+2,180-bend1,bend1_tol,constraint_type) ; cst_st += st
				
	# level 2
	if bend2:
		if nres >=10 and ( bulgepos_fromN >=5 and bulgepos_fromC <=-6 ):
			st = CstTypeAngleConstraints(pos1-4,pos1,pos1+5,180-bend2,bend2_tol,constraint_type) ; cst_st += st
			st = CstTypeAngleConstraints(pos2-5,pos2,pos2+4,180-bend2,bend2_tol,constraint_type) ; cst_st += st

	return cst_st



def HairpinPairingResidues(blue,segment1,segment2,str1_offset,str2_offset,bpos):
	# segment1 is the first strand in the hairpin according to topology
	s1 = blue.segment_dict[segment1]
	s2 = blue.segment_dict[segment2]
	map1={} ; map2={}
	pairs=[]

	if bpos == None:
		for i in range(len(s2.bp_data)):
			pos2 = int(s2.bp_data[i+str2_offset][0]) # Antiparallel
			if i < len(s1.bp_data) and s2.bp_data[i+str2_offset+1][2] != 'EA':
				pos1 = int(s1.bp_data[-1-i+str1_offset][0])
			else:
				break
			pairs.append([pos1,pos2])

	if bpos != None and Bulged(s1): 
		b1pos = bpos
		shift = int(s1.bp_data[-1][0]) - b1pos
		for i in range(len(s2.bp_data)-1):
			pos2 = int(s2.bp_data[i+str2_offset][0]) # Antiparallel
			if i < shift:
				pos1 = int(s1.bp_data[-1-i+str1_offset][0])
			elif i+1 < len(s1.bp_data):
				pos1 = int(s1.bp_data[-1-i-1+str1_offset][0])
			else:
				break
			pairs.append([pos1,pos2])
	# The bulgepos+1 is the one paired to the second strand

	if bpos != None and Bulged(s2):
		b2pos = bpos
		shift = b2pos - int(s2.bp_data[0][0])
		for i in range(len(s1.bp_data)):
			pos1 = int(s1.bp_data[-1-i+str1_offset][0])
			if i < shift:
				pos2 = int(s2.bp_data[i+str2_offset][0])
			elif i+1 < len(s2.bp_data):
				pos2 = int(s2.bp_data[i+1+str2_offset][0])
			else:
				break
		pairs.append([pos1,pos2])

	return pairs

def HbondsBulgedStrand(**kwargs):
	strand1 = kwargs.get('strand1')
	strand2 = kwargs.get('strand2')
	blue = kwargs.get('blueprint')
	str1_offset = kwargs.get('str1_offset', 0)
	str2_offset = kwargs.get('str2_offset', 0)
	b1pos = kwargs.get('bulge_position')

	pair1 = HairpinPairingResidues(blue,strand1,strand2,str1_offset,str2_offset,b1pos)
	seg1 = blue.segment_dict[strand1]
	seg2 = blue.segment_dict[strand2]

	if Bulged(seg1):
		seg = seg1
	elif Bulged(seg2):
		seg = seg2
	else:
		print('Warning: None of the strands is bulged')
#
#	b1pos = Bulged(seg)
	hblist=[]
	# to the left of the bulge
	i=0
	pos = b1pos
	while pos > seg.bp_data[0][0]:
		pos = (b1pos-2)-2*i
		i+=1
		hblist.append(pos)
    
	pos = b1pos
	# to the right of the bulge
	i=0
	while pos < seg.bp_data[-1][0]:
		pos = (b1pos+1)+2*i
		i+=1
		hblist.append(pos)
	hbpairs=[]    
	for pair in pair1:
		hbpair = list( set(hblist).intersection(set(pair)) )
		if len(hbpair)==1:
			hbpairs.append(pair)

	# get the type of beta-turn involved
	turnABEGO = ""
	for k in range(seg1.bp_data[-1][0], seg2.bp_data[0][0]-1):
		turnABEGO += blue.bp_data[k][2][1]
#	print turnABEGO
	cst_st = ''
	first = True
	for hbpair in hbpairs:
		pos1,pos2 = hbpair
		st = CircularHBondConstraints(pos1,pos2) ; cst_st += st
		if first == True:
			first = False
			if len(turnABEGO) == 2:
				st = CircularHBondConstraints(pos2,pos1) ; cst_st += st
			elif len(turnABEGO) == 3 and turnABEGO == "AAG":
				st = CircularHBondConstraints(pos2-1,pos1) ; cst_st += st
			elif len(turnABEGO) == 4 and turnABEGO == "AAAG":
				st = CircularHBondConstraints(pos2,pos1) ; cst_st += st
				st = CircularHBondConstraints(pos2-1,pos1) ; cst_st += st
				st = CircularHBondConstraints(pos2-2,pos1) ; cst_st += st
			elif len(turnABEGO) == 5 and turnABEGO == "AAAGG":
				st = CircularHBondConstraints(pos2,pos1) ; cst_st += st
				st = CircularHBondConstraints(pos2-1,pos1) ; cst_st += st
				st = CircularHBondConstraints(pos2-3,pos1) ; cst_st += st

		else:
			st = CircularHBondConstraints(pos2,pos1) ; cst_st += st
			# Add hbond for bulge position (which is not C-alpha paired)
	
		if (b1pos+1) in hbpair: # position b1pos+1 is in the pairing list (not the bulgepos, so we add the hbond here)
			# identify paired position
			paired_pos = list( set(hbpair).difference(set([b1pos+1])) )[0]
			st = CircularHBondConstraints(b1pos,paired_pos) ; cst_st += st

	if blue.bp_data[pos1 - 4][2] != "H":	
		cst_st += CircularHBondConstraints(pos2 + 2, pos1 -2)

	return cst_st

def HbondsBulgedPairWithOffset(**kwargs):
	strand1 = kwargs.get('strand1')
	strand2 = kwargs.get('strand2')
	blue = kwargs.get('blueprint')
	str1_offset = kwargs.get('str1_offset', 0)
	str2_offset = kwargs.get('str2_offset', 0)
	b1pos = kwargs.get('bulge_position')

	pair1 = HairpinPairingResidues(blue,strand1,strand2,str1_offset,str2_offset,b1pos)
	seg1 = blue.segment_dict[strand1]
	seg2 = blue.segment_dict[strand2]

	if Bulged(seg1):
		seg = seg1
	elif Bulged(seg2):
		seg = seg2
	else:
		print('Warning: None of the strands is bulged')
#
#       b1pos = Bulged(seg)
	hblist=[]
	# to the left of the bulge
	i=0
	pos = b1pos
	while pos > seg.bp_data[0][0]:
		pos = (b1pos-2)-2*i
		i+=1
		hblist.append(pos)
	i=0
	pos = b1pos
	# to the right of the bulge
	while pos < seg.bp_data[-1][0]:
		pos = (b1pos+1)+2*i
		i+=1
		hblist.append(pos)
	hbpairs=[]
	for pair in pair1:
		hbpair = list( set(hblist).intersection(set(pair)) )
		if len(hbpair)==1:
			hbpairs.append(pair)

	cst_st = ''
	for hbpair in hbpairs:
		pos1,pos2 = hbpair
		st = CircularHBondConstraints(pos1,pos2) ; cst_st += st
		st = CircularHBondConstraints(pos2,pos1) ; cst_st += st

		if (b1pos+1) in hbpair: # position b1pos+1 is in the pairing list (not the bulgepos, so we add the hbond here)
 			# identify paired position
			paired_pos = list( set(hbpair).difference(set([b1pos+1])) )[0]
			st = CircularHBondConstraints(b1pos,paired_pos) ; cst_st += st

	if blue.bp_data[pos1 - 4][2] != "H":
		cst_st += CircularHBondConstraints(pos2 + 2, pos1 -2)

	return cst_st

def AllSheetSegmentPairs(blue):
	cst_st=''
	pair1 = HairpinPairingResidues(blue,'E1','E2')
	pair2 = HairpinPairingResidues(blue,'E2','E3')
	pair3 = HairpinPairingResidues(blue,'E3','E4')

	pairs=[]
	pairs.extend(pair1)
	pairs.extend(pair2)
	pairs.extend(pair3)		

	dic_pairs={}
	for pair in pairs:
		a,b = pair
		seg_a = blue.residue_segment(a)
		seg_b = blue.residue_segment(b)
		dic_pairs.setdefault(a,{})
		dic_pairs.setdefault(b,{})
		dic_pairs[a][seg_b]=b
		dic_pairs[b][seg_a]=a
	
	return dic_pairs		

def CoordinateConstraints(**kwargs):
	blue = kwargs.get('blueprint')
	blue_positions = kwargs.get('blue_positions')
	pdb = kwargs.get('pdb')
	pdb_positions = kwargs.get('pdb_positions')
	chain = kwargs.get('chain')
	sd = kwargs.get('sd',0.5)
	d0=0
	# read pdb
	ca=[]
	for line in open(pdb):
		if line.split()[0] == 'ATOM' and line.split()[2]=='CA' and line.split()[4]==chain:
			coord = map(float,line.split()[6:9])         
			ca.append(coord)		
	st=''
	for i,j in zip(blue_positions,pdb_positions):
		x,y,z = ca[j]
		st += "CoordinateConstraint CA %i CA %i %.3f %.3f %.3f HARMONIC %.2f %.2f\n" %(i,i,x,y,z,d0,sd)

	return st
	
def HbondsRegularHairpin(**kwargs):
	strand1 = kwargs.get('strand1')
	strand2 = kwargs.get('strand2')
	blue = kwargs.get('blueprint')
	str1_offset = kwargs.get('str1_offset', 0)
	str2_offset = kwargs.get('str2_offset', 0)
	
	pair1 = HairpinPairingResidues(blue,strand1,strand2, str1_offset,str2_offset,None)
	seg1 = blue.segment_dict[strand1]
	seg2 = blue.segment_dict[strand2]

	inipos = seg1.bp_data[-1][0] # start counting from hairpin loop
	pos=inipos
	hblist=[pos]
	# all residues of seg1 are hbonded to seg2
	while pos >= seg1.bp_data[0][0]+2:
		pos -= 2
		hblist.append(pos)

	hbpairs=[]
	for pair in pair1:
		hbpair = list( set(hblist).intersection(set(pair)) )
		if len(hbpair)==1:
			hbpairs.append(pair)
	# get the type of beta-turn involved
	turnABEGO = ""
	for k in range(seg1.bp_data[-1][0], seg2.bp_data[0][0]-1):
		turnABEGO += blue.bp_data[k][2][1]

	first = True
	cst_st = ''
	for hbpair in hbpairs:
		pos1,pos2 = hbpair
		st = CircularHBondConstraints(pos1,pos2) ; cst_st += st
		if first == True:
			first = False
			if len(turnABEGO) == 2:
				st = CircularHBondConstraints(pos2,pos1) ; cst_st += st
			elif len(turnABEGO) == 3 and turnABEGO == "AAG":
				st = CircularHBondConstraints(pos2-1,pos1) ; cst_st += st
			elif len(turnABEGO) == 4 and turnABEGO == "AAAG":
				st = CircularHBondConstraints(pos2,pos1) ; cst_st += st
				st = CircularHBondConstraints(pos2-1,pos1) ; cst_st += st
				st = CircularHBondConstraints(pos2-2,pos1) ; cst_st += st
			elif len(turnABEGO) == 5 and turnABEGO == "AAAGG":
				st = CircularHBondConstraints(pos2,pos1) ; cst_st += st
				st = CircularHBondConstraints(pos2-1,pos1) ; cst_st += st
				st = CircularHBondConstraints(pos2-3,pos1) ; cst_st += st

		else:
			st = CircularHBondConstraints(pos2,pos1) ; cst_st += st

	# account for the last hbond between strand and a loop
	# excluding the tryptophan corner.
	if blue.bp_data[pos1 - 4][2] != "H":
		cst_st += CircularHBondConstraints(pos2 + 2, pos1 -2)
            
	return cst_st


def CutPDB(pdb,cut1,cut2,outpdb):	
	infile = pdb
	fileout = open(outpdb,'w')

	# cut in the range[cut1,cut2], both included
	at_index=0
	for line in open(infile):
        	if 'ATOM' in line:
                	res_index = int(line.split()[5])
	                if res_index < cut1:
        	                continue
	                elif res_index >= cut1 and res_index <= cut2:
        	                at_index+=1
                	        shift = cut1-1
                        	res_index -= shift
	                        line2 = line[:4] + '%7i' %(at_index) + line[11:22] + '%4i' %(res_index) + line[26:]
	                        fileout.write(line2)
        	        elif res_index > cut2:
                	        continue
	fileout.write('END')

