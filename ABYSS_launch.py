#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# @package ABYSS_launch.py
# @author Florian Charriat

"""
	The ABYSS_launch script
	=======================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 9/03/2018
	:version: 0.1

	Script description
	------------------

	This program is used to assemble all fasta files in a directory using the ABYSS tool. The assembly
	 will be done with different lengths of kmère (20, 30, 40, 50, 60, 70, 80 and 90). A new pipeline
	 has been created, please watch the Assembly_pipeline at https://github.com/FlorianCHA/AssemblyAndAnnotation_pipeline

	Example
	-------

	>>> ABYSS_launch.py -d /homedir/user/work/data -o /homedir/user/work/result

	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						Display ABYSS_launch.py version number and exit
						
	Input mandatory infos for running:
		- \-d <path/to/directory>, --directory <path/to/directory>
						path of directory that contains all the fastq files which must be assembled
						
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						path of the output directory

"""


########## Module ###############
## Python modules
import argparse, os, sys

#Import module_Flo
from module_Flo import verifDir, createDir , form, header_script, footer_script



if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to assemble all fastq files in a directory using \
	the ABYSS tool. The assembly will be done with different lengths of kmère (20, 30, 40, 50, 60, 70, 80 and 90) ''')
	parser.add_argument('-v', '--version', action='version', version=f'You are using ABYSS_launch version: {version}', help='Display ABYSS_launch version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-d', '--directory',type = str, required=True, dest = 'dirPath', help = 'Path of directory that contains all the fasta files which must be assembled')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	directory = os.path.abspath(args.dirPath)
	outDir= os.path.abspath( args.outdirPath)


########### Gestion directory ##############
	directory = verifDir(directory,True)
	outDir = verifDir(outDir)
	bash = f'{outDir}script_bash'
	name_directory = [outDir,f'{outDir}error_files', f'{outDir}out_files',bash,f'{outDir}result']
	createDir(name_directory)

############### start message ########################
	header_script('ABYSS_launch',version)


############# Main #########################
	nbJob = 0
	nbGenome = 0
	# Loop on all file in the input directory which contains all fastq
	for file in os.listdir(directory):
		# Filters files with the correct extension
		if file.endswith('_R1.fastq.gz')==True :
			nbGenome += 1
			isolate=file.replace('_R1.fastq.gz','')
			print(form('\nLancement des jobs pour : '+isolate+'\n','green',['bold','underline']))
			# Loop on kmere parameters list for ABySS (This pipeline launch all assembly with a kmere of 20,30,40,50,60,70,80 and 90.
			for kmers in [20,30,40,50,60,70,80,90]:
				# Write a sh script for AGAP Cluster with log output, log error etc ....
				nameScript =bash+'/abyss_assembly_'+isolate+'_'+str(kmers)+'.sh'
				SCRIPT=open(nameScript,'w')
				SCRIPT.write('#$ -o '+outDir+'out_files/abyss_assembly_'+isolate+'_'+str(kmers)+'.out\n#$ -e '+outDir+'error_files/abyss_assembly_'+isolate+'_'+str(kmers)+'.err\nmodule load bioinfo/abyss/1.9.0;\n')
				SCRIPT.write('echo $PATH\n')
				SCRIPT.write('mkdir -p '+outDir+'result/'+isolate+'/abyss_assembly_'+isolate+'_'+str(kmers)+';\n')
				SCRIPT.write('cd '+outDir+'result/'+isolate+'/abyss_assembly_'+isolate+'_'+str(kmers)+';\n')
				SCRIPT.write("/usr/local/bioinfo/abyss/1.9.0/bin/abyss-pe name="+isolate+"_"+str(kmers)+" k="+str(kmers)+" in='"+directory+file+" "+directory+file.replace('_R1','_R2')+"' -o abyss_assembly_"+isolate+"_"+str(kmers)+".fasta;\n")
				SCRIPT.close()
				# Launch all script
				os.system('qsub -l mem_free=50G -l h_vmem=60G -q normal.q -V -cwd '+nameScript)
				nbJob += 1


############## end message ###########################

	footer_script()



