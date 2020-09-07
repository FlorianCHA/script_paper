#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package secretome_Pipeline.py
# @author Florian Charriat

"""
	The secretome_Pipeline script
	=============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 18/04/2018
	:version: 0.1

	Script description
	------------------

	This program is used to predict secretome with SignalP, TagetP, Phobius.\n
	This program uses the comparaisonSecretome script to retrieve and compare information from the secretome prediction tools.\n
	And the selection_TMHMM script to select, from TMHMM ouput, only protein with no TH domain or only one TH in 60 first aa.\n
	This program uses also the eliminateREmotif script to eliminate, from ps_scan ouput, the protein with a RE retention motif.\n
	Please make sure that this script is present in the directory.
	
	Example
	-------

	>>> secretome_Pipeline.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display secretome_Pipeline.py version number and exit
						
	Input mandatory infos for running:
		- \-f <path/to/fasta/file>, --file <path/to/fasta/file>
						path of fasta files that contains all the protien of strain
		- \-p <path/to/prosite.dat/file>, --outdirPath <path/to/prosite.dat/file>
						Path of prosite.dat file. You can upload the file at ftp://ftp.expasy.org/databases/prosite/prosite.dat						
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						Path of the output directory
		- \-fo                             --force
						force the script to remove the old output data
"""


########## Module ###############
## Python modules
import argparse, os, sys,header_script

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId, verifFichier, header_script,footer_script

#Prosite file 
#PathPrositeFile = '/homedir/charriat/BioInfo_Tools/ps_scan/prosite.dat'

if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to predict secretome. This program uses :

	- SignalP, TagetP, Phobius for detect the presence of signal peptide cleavage sites.
	- ComparaisonSecretome script to retrieve and compare information from the secretome prediction tools.
	- Selection_TMHMM script to select, from TMHMM ouput, only protein with no TH domain or only one TH in 60 first aa.
	- EliminateREmotif script to eliminate, from ps_scan ouput, the protein with a RE retention motif.
	
Please make sure that this script is present in the directory.''',formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display secretome_Pipeline version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-f', '--file',type = str, required=True, dest = 'dirPath', help = 'path of fasta files that contains all the protein of strain')
	filesreq.add_argument('-p', '--prosite',type = str, required=True, dest = 'prositePath', help = 'Path of prosite.dat file. You can upload the file at ftp://ftp.expasy.org/databases/prosite/prosite.dat')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')

	filesNoreq = parser.add_argument_group('Input not mandatory infos for running')
	filesNoreq.add_argument('-fo', '--force', action='store_true', dest = 'force', help = 'force the script to remove output data')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	fasta_file = os.path.abspath(args.dirPath)
	outDir= os.path.abspath(args.outdirPath)
	force = args.force
	prosite_dat =  os.path.abspath(args.prositePath)
	verifFichier(prosite_dat)
########### Gestion directory ##############
	verifFichier(fasta_file)
	outDir = verifDir(outDir)
	bash = outDir+'script_bash'
	name_directory = [outDir,outDir+'error_files', outDir+'out_files',bash]
	createDir(name_directory)


############### start message ########################
	header_script('secretome_Pipeline',version)
############# Main #########################
	nbfile = 0
	file = fasta_file.split('/')[-1]
	if isFasta(file):
		idFile = recupId(file)
		nbfile +=1
		fasta = open(fasta_file,'r')
		lines = fasta.readlines()
		fasta.close()
		nb = 0 # Initiate a variable for split fasta into multi fasta with less than 400 sequences
		nbSeq = 0 # Initiate a variable for know the number of sequence present in the fasta
		part = 1 # Initiate a variable for rename files to differentiate fasta files split

		########## Create output directory ############
		outDir_result = f'{outDir}Result/{idFile}'
		if os.path.exists(outDir_result) == True and force == False:
			raise ValueError(form(f"Output directory: '{outDir_result}' already exist, please remove it or use the -- force option for remove automatically the directory","red","bold"))
		if os.path.exists(outDir_result) != 0 and force == True :
			os.system(f'rm -r {outDir}Result/{idFile}')
		outDirFasta = f'{outDir_result}/0_fasta-files'
		outDir_comparaison = f'{outDir_result}/1_predicted/'
		outDir_selectTMHMM = f'{outDir_result}/2_TMHMM/'
		outDir_PS_scan = f'{outDir_result}/3_PS-scan/'
		outDir_data_final = f'{outDir_result}/4_data-final/'
		outDir_effectorP = f'{outDir_result}/5_effectorP/'
		name_directory = [f'{outDir}Result',outDirFasta,outDir_comparaison,outDir_selectTMHMM,outDir_PS_scan,outDir_data_final,outDir_effectorP]
		createDir(name_directory)
			
		######### Create results files ##########
		outTargetP = f'{outDir_comparaison}{idFile}_targetP.txt'
		outputPhobius = f'{outDir_comparaison}{idFile}_phobius.txt'
		outputSignalP = f'{outDir_comparaison}{idFile}_signalP.txt'
		os.system(f'touch {outTargetP} {outputPhobius} {outputSignalP}')

		listeTargetp = []
		listePhobius = []
		listeSignalP = []
	
		for line in lines:
			# Split the fasta into multi fasta with max 400 sequences (format input for tools)
			if nb == 400 and line[0] == '>' :
				listePhobius.append(f'phobius.pl -short {outDirFasta}/{idFile}_part{part}.fasta  >> {outputPhobius};\n')
				listeTargetp.append(f'targetp -N {outDirFasta}/{idFile}_part{part}.fasta >> {outTargetP};\n')
				listeSignalP.append(f'signalp -u 0.24 -U 0.24 {outDirFasta}/{idFile}_part{part}.fasta >> {outputSignalP};\n')
				nb = 1
				nbSeq +=1
				part +=1
				f = open(f'{outDirFasta}/{idFile}_part{part}.fasta','a')
				f.write(line)
				f.close()
					
				
			elif line[0] == '>' :
				nb +=1
				nbSeq +=1
				f = open(f'{outDirFasta}/{idFile}_part{part}.fasta','a')
				f.write(line)
				f.close()
				
			else : 
				f = open(f'{outDirFasta}/{idFile}_part{part}.fasta','a')
				f.write(line)
				f.close()

		#Add all command in a script file for each fasta
		if nbSeq%400 != 0 :
			listePhobius.append(f'phobius.pl -short {outDirFasta}/{idFile}_part{part}.fasta  >> {outputPhobius};\n')
			listeTargetp.append(f'targetp -N {outDirFasta}/{idFile}_part{part}.fasta >> {outTargetP};\n')
			listeSignalP.append(f'signalp -u 0.34 -U 0.34 {outDirFasta}/{idFile}_part{part}.fasta >> {outputSignalP};\n')
					
			
		with open(f'{bash}/{idFile}_secretomeTools.sh','w') as f :

			f.write('\n########### Launch SignalP ###################\n\n')
			for elt in listeSignalP :
				f.write(elt)
			f.write('\n########### Launch targetP ###################\n\n')
			for elt in listeTargetp :
				f.write(elt)
			f.write('\n########### Launch Phobius ###################\n\n')
			for elt in listePhobius :
				f.write(elt)
		with open open(f'{bash}/{idFile}.sh','w') as f :
			f.write(f'#$ -e {outDir}error_files\n#$ -o {outDir}out_files\n#$ -l mem_free=10G\n#$ -N P_{idFile}_secretome\n#$ -q normal.q\n#$ -V\n\nmodule load bioinfo/signalp/4.1\nmodule load system/java/jre8\n\n')
			f.write(f'\n\n######### Launch Phobius, TargetP and SignalP #########\n\n')
			f.write(f'bash {bash}/{idFile}_secretomeTools.sh\n\n')
			f.write(f'######### Comparaison des 3 outils de prédiction #########\n\n')
			f.write(f'comparaisonSecretome.py -o {outDir_comparaison} --phobius {outputPhobius} --targetp {outTargetP} --signalp {outputSignalP} --rank 2 --fasta {fasta_file}\n')

	################################### TMHMM ##########################################
			f.write('\n\n######### Launch TMHMM #########\n\n')
			f.write(f'tmhmm -short {outDir_comparaison}{idFile}_secreted_1.fasta > {outDir_selectTMHMM}{idFile}_TMHMM.txt\n')
			f.write('\n\n######### Filtering proteins according to the TMHMM #########\n\n')
			f.write(f'selection_TMHMM.py -t {outDir_selectTMHMM}{idFile}_TMHMM.txt -f {outDir_comparaison}{idFile}_secreted_1.fasta -o {outDir_selectTMHMM}\n')
			
############################ PS_scan for RE retention motif ##########################
			f.write('\n\n######### Filtering proteins according to the pattern of retention in the ER with PS-scan #########\n\n')
			f.write(f'ps_scan.pl -o pff -p PS00014 -d {prosite_dat} {outDir_selectTMHMM}{idFile}_secreted_2.fasta > {outDir_PS_scan}{idFile}_ps_scan.txt\n')
			f.write(f'elimateREmotif.py -p {outDir_PS_scan}{idFile}_ps_scan.txt -f {outDir_selectTMHMM}{idFile}_secreted_2.fasta -o {outDir_PS_scan}\n')
			f.write(f'cp {outDir_PS_scan}{idFile}_secreted_3.fasta {outDir_data_final}{idFile}_secreted.fasta\n')
			print(form(f'Script create for {idFile}\n','green','bold'))

		
########################### Launch EffectoP #############################
			f.write('EffectorP.py -o {outDir_effectorP}{idFile}_stat.txt -E {outDir_effectorP}{idFile}.fasta -i {outDir_data_final}{idFile}_secreted.fasta')

############## summary message #######################

	print(form('\n-----------------------------------------------------------------------------------------------------------------------','red','bold'))
	print(form('Execution summary:\n','green',['bold','underline']))
	print(f'\n\tInput : \n\t\t- {fasta_file}')
	print('\n\tOutput :')
	print(f'\t\t - Résultat des prédictions des secretomes : {outDir_result}')
	print(f'\nThe secretome_Pipeline processed the strain {idFile}, please type the following command to run the created scripts :\n\n\t\t\t\t'+form(f'qsub {bash}/{idFile}.sh\n','green','bold'))
	print(form('----------------------------------------------------------------------------------------------------------------------','red','bold'))





############## end message ###########################
	footer_script()




