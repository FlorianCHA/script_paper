#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package alnToSnp.py
# @author Florian CHARRIAT



##################################################
## Modules
##################################################
#Import 
import sys, os
# Import module_Flo
from module_Flo import verifDir, createDir, form, verifFichier, fasta2dict, openfile, sort_human,header_script,footer_script
from packages.progress.bar import ChargingBar

## Python modules
import argparse, collections

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
##################################################
## Functions


##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	version = '0.1'

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme is used to keep only SNP from a alignement''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')

	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-a', '--alignement', type = str, required=True, dest = 'file', help = 'path of alignement file')
	filesreq.add_argument('-o', '--output',type = str, required=True, dest = 'output', help = 'path of output file')
	filesreq.add_argument('-t', '--type', type=str, required=True, dest = 'type', help = 'Alignement type')

######### Retrieve arguments ###########
	args = parser.parse_args()
	file = os.path.abspath(args.file)
	output = os.path.abspath(args.output)
	type =  args.type

########### Directory control ##############
	verifFichier(file)

############### start message ########################
	header_script('alnToSnp',version)
########### Main #####################################
	print("\nRetrieving the alignment file\n")
	# Use Biopython module for parse alignement
	alignment = AlignIO.read(file, type)
	nbSeq = len(alignment)
	length_Seq = len(alignment[1,:].seq)
	bar = ChargingBar("In process", max=length_Seq, suffix='%(percent)d%%')
	nb = nbG = nbN = nbGap = 0
	dico_count = collections.defaultdict(int)
	for i in range(0,length_Seq) :
		# Process position by position
		str_columns = alignment[:,i].upper()
		nb += 1
		# Eliminate position with at least one N in the sequence
		if 'N' in str_columns :
			nbN += 1
		# Eliminate position with at least one ga in the sequence
		elif '-' in str_columns :
			nbGap += 1
		# Checks if this position corresponding to a SNP
		if nbSeq != str_columns.count(str_columns[0]) and 'N' not in str_columns and '-' not in str_columns:
			nbG += 1
			# If it's a SNP, add this position to the new alignement object
			try :
				edited += alignment[:, i:i+1]
			except NameError:
				edited = alignment[:, i:i+1]
		elif nbSeq == str_columns.count(str_columns[0]) :
			dico_count[str_columns[0]] += 1
		bar.next()
	bar.finish()
	# Write the new alignement
	AlignIO.write(edited,output, 'fasta')
	print(form(f'The script kept {nbG} nucleotide of the {nb} total nucléotide of the alignement (Nb Gap = {nbGap}, Nb N = {nbN}, Nb A = {dico_count["A"]}, Nb C = {dico_count["C"]}, Nb G = {dico_count["G"]}, Nb T = {dico_count["T"]} )','green','bold'))
	with open(f'{output.split(".")[0]}_stat.out','w') as stat_file :
		stat_file.write(f'The script kept {nbG} nucleotide of the {nb} total nucléotide of the alignement (Nb Gap = {nbGap}, Nb N = {nbN}, Nb A = {dico_count["A"]}, Nb C = {dico_count["C"]}, Nb G = {dico_count["G"]}, Nb T = {dico_count["T"]} )')

############## end message ###########################
	footer_script()
