#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package Correction_blast.py
# @author Florian Charriat

########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human,indexEgale,indexDif,functionSens,comparaisonListe,isIn,selectBlast,header_script,footer_script

if __name__ == "__main__":

	version="0.1"

############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used  to correct orthologue groups of a sequence from fasta sequence and orthofinder output''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display Correction_blast version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'directory', help = 'Path of the directory which contains all OG fasta file')
	filesreq.add_argument('-c', '--count',type = str, required=True, dest = 'count', help = 'Path of the count result of orthofinder (format csv)')
	filesreq.add_argument('-db', '--database', type=str, required=True, dest='db',
						  help='Path of the all assembly merge file with scafold name : Souche_Scafold_1')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdir', help = 'Path of the output directory')



######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.directory)
	db = os.path.abspath(args.db)
	counts = os.path.abspath(args.count)
	outdir =  os.path.abspath(args.outdir)
########### Gestion directory ##############
	verifFichier(counts)
	directory = verifDir(directory, True)
	outdir = verifDir(outdir)
	# Create all directory output
	createDir([outdir,outdir+'OG_fasta',outdir+'blast_result',outdir+'script'])


############### start message ########################
	header_script('Correction_blast',version)


############# Main ##################################

	listeOG = []
	i = 1
	# Parse fasta OG one by one
	for elt in os.listdir(directory) :
		# Checks if it's a fasta file
		if elt.endswith('.fasta') :
			print(form(f'Traitement du fichier {elt} : ', 'green', 'bold'))
			print()
			listeOG.append(elt.split(".")[0])
			print(form(f'\t- Open file', 'white', 'bold'))
			OG_file = fasta2dict(f'{directory}/{elt}')
			length_prev = 0
			# This part select the smallest sequence of the fasta for blast
			for name in OG_file.keys():
				length = len(OG_file[name].seq)
				# This condition select the smallest sequence of the fasta for blast and remove of the selection, the outgroup of the data set
				# and sequence with N
				if length_prev > length or length_prev == 0 and 'XX' not in str(OG_file[name].seq):
					if 'BR29' not in name and 'Pm1' not in name and 'DsLIZ' not in name and 'MGG' not in name and 'XX' not in str(
							OG_file[name].seq):
						min = [name, OG_file[name].seq]
						length_prev = length
			print(form(f'\t- Lancement du Blast', 'white', 'bold'))
			# Write the query for the blast analysis with the smallest sequence
			with open(f'{outdir}/OG_fasta/{elt}', 'w') as og_file:
				record = SeqRecord(min[1], id=elt.split(".")[0], name=elt.split(".")[0],
								   description=f'| length : {len(min[1])}')
				SeqIO.write(record, og_file, "fasta")
			# The time to make blast analysis one by one is long, if you want earn time, you can comment the os.system line and all other parts (Retrieve strain, parse blast and write output)
			os.system(f'blastn -query {outdir}OG_fasta/{elt} -db {db} -evalue 1e-4 > {outdir}blast_result/{elt.split(".")[0]}_blast.txt')
			# Then decomment this two line for create all sh script. Launch all script in parallel then remove comment only of 3 next part !
			# with open(f'{outdir}script/launcher_{i}.sh','w') as f :
			# 	f.write(f'#!/bin/bash\nmodule load bioinfo/ncbi-blast/2.6.0\nblastn -query {outdir}OG_fasta/{elt} -db {db} -evalue 1e-4 > {outdir}blast_result/{elt.split(".")[0]}_blast.txt')  # "6 qseqid sseqid pident qlen length mismatch gapope evalue qstart qend  sstart send "
			i += 1
			print()
			print(form('-' * 100, 'white', 'bold'))
			print()
	# ################################## Retrieve strain absent of OG #########################
	nb = 0
	dico = {}
	print("Retrieve strain absent of OG")
	with open(counts,'r') as file_count,\
			open(f'{outdir}gene.count_prov.csv','w') as result :
		header = f'Groupe{file_count.readline()}'.replace('\tTotal','')
		listeProtein = header.replace('\n','').split('\t')[1:len(header.split('\t'))]
		result.write(header)
		for line in file_count :
			name = line.split()[0]
			if name in listeOG :
				nb += 1
				total = line.split()[-1]
				new_line = line.replace(f'\t{total}', '')
				result.write(new_line)
				OG = line.split('\t')[0]
				count = line.split('\t')[1:-1]
				if '0' not in count :
					print(f'{OG} is in core genome')
				else :
					listeIndex = indexEgale(count, '0')
					listeCount = []
					print(OG)
					for i in listeIndex:
						listeCount.append(listeProtein[i])
					dico[OG] = listeCount

	################################## Parse blast result #######################################
	dico_blast = {}
	print("Parse blast result")
	for name in os.listdir(f'{outdir}blast_result/'):
		print(name)
		OG = name.split('_')[0]
		if OG in dico.keys() :
			# Use a function which select only hit with 90 percent of identity and coverage as well as a max intron of 1000 pb
			dico_blast[OG] =  [elt.split('_')[0] for elt in selectBlast(f'{outdir}blast_result/{name}',90,1000) if elt.split('_')[0] in dico[OG]]

	############################ Write new OG_count output ############################################
	print("Write new OG_count output")
	with open(counts,'r') as count_file ,open(f'{outdir}gene.count.csv','w') as result:
		header = count_file.readline()
		listeProtein = header.replace('\n', '').replace('Groupe','').split('\t')
		result.write(header)
		for line in count_file:
			OG = line.split('\t')[0]
			lineSplit = line.replace('\n','').split('\t')
			print('Before : ',lineSplit)
			if OG in dico_blast.keys():
				print(OG)
				for elt in dico_blast[OG]:
					# This function retrieve the position in the line of the element to change
					listeI = indexEgale(listeProtein, elt)
					for i in listeI:
						print(i)
						print(len(lineSplit))
						# Add a 999 score for all presence by blast correction
						lineSplit[i] = '999'
			print("After : ",lineSplit)
			newline = '\t'.join(lineSplit)
			result.write(newline+'\n')

	############## end message ###########################
	footer_script()




