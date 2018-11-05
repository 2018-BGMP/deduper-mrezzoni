#!/usr/bin/env python3


# module load python3/3.6.1
import argparse
import re

def get_arguments():
	parser = argparse.ArgumentParser(description="Allows users to know how to navigate this program")
	parser.add_argument("-f", "--file", help="Absolute file path. Must include .sam in the name", required=True, type=str)
	parser.add_argument("-p", "--paired", help="(Currently unsupported) Designates file is paired end (not single-end)", required=False, type=str)
	parser.add_argument("-u", "--umi", help="Designates file containing the list of UMIs", required=True, type=str)
	return parser.parse_args()
	
################################
###### Argparse Variables ######
################################

args = get_arguments()
file = args.file
paired = args.paired
umi = args.umi

if paired != None:
    print("I suck and did not take the time to make this script capable of accepting paired end reads.")
    sys.exit()

#######################################
######## High Level Function ##########
#######################################
	
def correct_5prime_start(FLAG:str, POS:str, CIGAR:str) -> int:
	'''
	FLAG: Bitwise flag indicative of forward/reverse strand
	POS: Left-most mapping position
	CIGAR: Cigar string
	
	returns: Corrected left-most mapping position for forward (sense) strand or
	right-most mapping position for reverse (antisense) strand. Returns original 
	position if there are no cigar operations.
	'''
	if((int(FLAG) & 16) != 16):
		# search the forward (sense) strand
		sense_S = re.findall("^\d+S", CIGAR) # search beginning of string for numbers that preceed S
		if sense_S:
			# grab any softclipping isolated by findall
			sense_clipped = int(sense_S[0][:-1]) # isolate the number
			new_pos = int(POS) - sense_clipped # adjust the position for sense strand
			return(new_pos)
		else: 
			return(int(POS))
		
	else:
		# search the reverse (antisense) strand
		new_pos = int(POS)
		anti_S = re.findall("\d+S$", CIGAR) # search end of string for numbers that preceed S
		if anti_S != []:
			sum_S = 0
			S_list = []
			for soft_clip in anti_S:
				# grab any softclipping isolated by findall
				soft_clips = soft_clip[:-1] 
				soft_clips = int(soft_clips)
				S_list.append(soft_clips)
				sum_S = sum(S_list)
			new_pos += sum_S
		
		M = re.findall("\d+M", CIGAR)
		if M != []:
			sum_M = 0
			M_list = []
			for match in M:
				# grab any matches/mismatches that were isolated by findall
				matches = match[:-1]
				matches = int(matches)
				M_list.append(matches)
				sum_M = sum(M_list)
				new_pos += sum_M
		
		D = re.findall("\d+D", CIGAR)
		if D != []:
			sum_deletions = 0
			deletion_list = []
			for deletion in D:
				# grab any deletions that were isolated by findall
				deletions = deletion[:-1]
				deletions = int(deletions)
				deletion_list.append(deletions)
				summed_deletions = sum(deletion_list)
				new_pos += sum_deletions
				
		N = re.findall("\d+N", CIGAR)
		if N != []:
			sum_skips = 0
			skip_list = []
			for skip in N:
				# grab any skipped regions due to alternative splicing that were isolated by findall
				skips = skip[:-1]
				skips = int(skips)
				skip_list.append(skips)
				summed_skips = sum(skip_list)
				new_pos += summed_skips
				
		return(new_pos)

#######################################
########## Open Output File ###########
#######################################	

output_name = re.search("([A-Za-z0-9]+\.sam)$", file).group(1) # removes the full path to preserve original .sam file name
deduped = open(output_name+"_deduped", "wt")

#######################################
######## Check for UMI list ###########
#######################################	

if umi != None:
	with open(umi, "r") as umis:
		umi_list = []
		for u in umis:
			umi_list.append(u.strip())
else:
	print("This code is currently only equipped to handle UMIs. Randomers will be a feature in the next update.")
	sys.exit()
	
#######################################
##### Open Input, Grab Uniques ########
#######################################	

uniq_reads = set()
chroms_seen = set()
total_ct = 0
dup_ct = 0

with open(file, "rt") as f:
	while True: # keep looping until the break

		full_sam = f.readline() # maintains entire SAM line
		sam_list = full_sam.strip().split() # splits into list to isolate components of interest
		
		if sam_list == []:
			break
		
		if full_sam.startswith("@"):
			# write out header lines
			deduped.write(full_sam)
			continue
		
		total_ct += 1
		
		umi = sam_list[0][-8:]
		strand = sam_list[1]
		chrom = sam_list[2]
		left_most = sam_list[3]
		cig_string = sam_list[5]
		
		if total_ct == 1:
			# special case for the very first the time loop runs
			chroms_seen.add(chrom)
		
		if umi in umi_list:

			new_start = correct_5prime_start(strand, left_most, cig_string) # make these variables
			#print(new_start)
			
			set_sam = (umi, strand, chrom, new_start) # umi, strand, chromosome, left most position
			#print(set_sam)
			
			if set_sam not in uniq_reads:
				# must be a unique read - immediately write out to new file
				uniq_reads.add(set_sam)
				deduped.write(full_sam)
			
			else:
				dup_ct +=1
				continue
			
			if chrom not in chroms_seen:
				# refresh the set to conserve memory once an entire chromosome has been seen
				chroms_seen.add(chrom)
				uniq_reads = set()
				
				deduped.write(full_sam) # first unique instant
				uniq_reads.add(set_sam) # add it to the set immediately
		
		else:
		# dont include UMIs that are not legitimate
			continue

deduped.close()

################################
## Relevant Output Statements ##
################################

dup_proportion = (dup_ct / total_ct) * 100
dup_proportion = round(dup_proportion, 2)
print("PCR duplicates compose {}% of your file".format(dup_proportion))



