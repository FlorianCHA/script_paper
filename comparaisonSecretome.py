#!/usr/local/bioinfo/python/3.4.3_build2/bin/python
# -*- coding: utf-8 -*-
# @package comparaisonSecretome.py
# @author Florian Charriat

"""
	The comparaisonSecretome script
	===============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 18/04/2018
	:version: 0.1

	Script description
	------------------

	This program is used to retrieve and compare information from the output of secretome prediction tools (signalP, targetP and Phobius). This program is used by secretome_Pipeline
	
	Example
	-------

	>>> comparaisonSecretome.py  --signalp /homedir/user/work/data/output_signalp.txt --targetp /homedir/user/work/data/output_targetp.txt --phobius /homedir/user/work/data/output_phobius.txt -o /homedir/user/work/result
	
	>>> comparaisonSecretome.py --signalp /homedir/user/work/data/output_signalp.txt -o /homedir/user/work/result


	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display comparaisonSecretome.py version number and exit
						
	Input mandatory infos for running:
								
		- \-s <path/to/output/signalp/file>, --signalp <path/to/output/signalp/file>
						path of the signalP output file
		- \-t <path/to/output/targetp/file>, --targetp <path/to/output/targetp/file>
						path of the targetp output file
		- \-p <path/to/output/phobius/file>, --phobius <path/to/output/phobius/file>
						path of the phobius output file		
		- \-f <path/to/output/phobius/file>, --fasta <path/to/fasta/file>
						path of the phobius output file		
		- \-r <int>, --rank <int>
						rank mini by default rank = 3, here rank 1 = protein secreted detect by all tools, rank2 = protein secreted detect by two tools, rank 3 : protein secreted detect by one tool.	
		- \-o <path/to/output/directory>, --outdirPath <path/to/output/directory>
						path of the output directory

"""


########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human


if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''This program is used to predict secretome with SignalP, TagetP and Phobius''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display comparaisonSecretome version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-s', '--signalp',type = str, default = 'None', dest = 'signalp', help = 'Path of the signalP output file')
	filesreq.add_argument('-t', '--targetp',type = str, default = 'None', dest = 'targetp', help = 'Path of the targetp output file')
	filesreq.add_argument('-p', '--phobius',type = str, default = 'None', dest = 'phobius', help = 'Path of the phobius output file')
	filesreq.add_argument('-f', '--fasta',type = str,  required=True, dest = 'fasta', help = 'path of the phobius output file')
	filesreq.add_argument('-r', '--rank',type = int, default =3, dest = 'rank', help = 'rank mini by default rank = 3, here rank 1 = protein secreted detect by all tools,\nrank2 = protein secreted detect by two tools,\n rank 3 : protein secreted detect by one tool.	')
	filesreq.add_argument('-o', '--outdir',type = str, required=True, dest = 'outdirPath', help = 'Path of the output directory')

	
######### Recuperation arguments ###########
	args = parser.parse_args()
	signalp = args.signalp
	targetp = args.targetp
	phobius = args.phobius
		
	if args.signalp != "None" :
		signalp = os.path.abspath(args.signalp)
		verifFichier(signalp)
		

	if args.targetp != "None" :
		targetp = os.path.abspath(args.targetp)
		verifFichier(targetp)

	if args.phobius != "None" :
		phobius = os.path.abspath(args.phobius)
		verifFichier(phobius)
		
	if signalp == "None" and targetp == "None" and phobius == "None":
		raise ValueError(form("ERROR : You did not give an output file from one or more secretarial prediction tools.\n Please use the --phobius | --targetp | --signalp parameters,\
		 for more information look at the help (parserSecretome -h)",'red','bold'))
		
	outDir= os.path.abspath(args.outdirPath)
	outDir = verifDir(outDir)
	fasta = os.path.abspath(args.fasta)
	Id = recupId(fasta.split('/')[-1])
	rankMini = args.rank

########### Gestion directory ##############

	name_directory = [outDir]
	createDir(name_directory)
	

############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("     Welcome in comparaisonSecretome (Version " + version + ")     ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

####################################### Main #################################################
	listesignalp = []
	listetargetp = []
	listephobius = []
	dico = {}
	listeFinal = []
	GeneSignalP = []
	GeneTargetP = []
	GenePhobius  = []
	
######################### Analyse signalp file ################
	if signalp != "None" :
		print(form(f"- Analyze the  SignalP output ({signalp}) .\n","green","bold"))
		file_signalp = open(signalp,'r')
		lines = file_signalp.readlines()
		file_signalp.close()
		
		for line in lines :
			if line[0] != '#' :
				line = line.split()
				if line[9] == 'Y':
					listesignalp.append(line[0])
					gene = line[0].split('.t')[0]
					if gene not in GeneSignalP:
						GeneSignalP.append(gene)

######################### Analyse targetP file ################

	if targetp != "None" :
		print(form(f"- Analyze the  TargetP output ({targetp}).\n","green","bold"))
		file_targetp = open(targetp,'r')
		lines = file_targetp.readlines()
		file_targetp.close()
		start = False # Initiate the variable which determined if "--------" seq is in start of end of the output table
		for line in lines :
			if '----------------------------' in line and debut == False :
				start = True
			elif '----------------------------' in line and debut == True :
				start = False
			elif start :
				line = line[:-1].split()
				if line[5] == 'S':
					listetargetp.append(line[0])
					gene = line[0].split('.t')[0]
					if gene not in GeneTargetP:
						GeneTargetP.append(gene)


	

######################### Analyse phobius file ################

	if phobius != "None" :
		verifFichier(phobius)
		print(form(f" - Analyze the  phobius output ({phobius}).\n","green","bold"))
		file_phobius = open(phobius,'r')
		lines = file_phobius.readlines()
		file_phobius.close()
		
		for line in lines :
			if 'SEQENCE ID' not in line : # For remove header of the phobius output
				line = line[:-1].split()
				dico[line[0]] = [4,'No','No','No']
				if line[2] == 'Y' :
					listephobius.append(line[0])
					gene = line[0].split('.t')[0]
					if gene not in GenePhobius:
						GenePhobius.append(gene)
					
######################## Parse all result into a dico ##########################""

	# Initialise les variables qui vont permet de connaire le nombre de protéines secrété prédite par les 3 outil ou seulement deux ou un seul
	# Initiate variable for give a 'rank' to the secreted protein to know if they were predicted by 1, 2 or 3 tools
	listeGene = []
	listeRank1 = []
	nbRank1 = 0
	listePhTar = [] 
	nbPhTar = 0
	listePhSign = []
	nbPhSign = 0 
	listeSignTar = []
	nbSignTar = 0 
	listeTar = [] 
	nbTar = 0 
	listePh = []
	nbPh = 0 
	listeSign = []
	nbSign = 0
	for elt in dico.keys() :
		gene = elt.split('.t')[0]
		if gene not in listeGene :
			listeGene.append(gene)
		if elt in listesignalp:
			dico[elt][1] = 'Yes'
			dico[elt][0] = dico[elt][0] - 1
			
		if elt in listetargetp:
			dico[elt][2] = 'Yes'
			dico[elt][0] = dico[elt][0] - 1
			
		if elt in listephobius:
			dico[elt][3] = 'Yes'
			dico[elt][0] = dico[elt][0] - 1
		if  dico[elt][1] == 'Yes' and dico[elt][2] == 'Yes' and dico[elt][3] == 'Yes' :
			if gene not in listeRank1 :
				listeRank1.append(gene)
			nbRank1 += 1
		elif  dico[elt][2] == 'Yes' and dico[elt][3] == 'Yes':
			if gene not in listePhTar :
				listePhTar.append(gene)
			nbPhTar += 2
		elif  dico[elt][1] == 'Yes' and dico[elt][3] == 'Yes':
			if gene not in listePhSign :
				listePhSign.append(gene)
			nbPhSign += 1
		elif  dico[elt][1] == 'Yes' and dico[elt][2] == 'Yes':
			if gene not in listeSignTar :
				listeSignTar.append(gene)
			nbSignTar += 1
		elif  dico[elt][1] == 'Yes' :
			if gene not in listeSign :
				listeSign.append(gene)
			nbSign += 1	
		elif  dico[elt][2] == 'Yes' :
			if gene not in listeTar :
				listeTar.append(gene)
			nbTar += 1
		elif  dico[elt][3] == 'Yes' :
			if gene not in listePh :
				listePh.append(gene)
			nbPh += 1

			
	nbRank1 = 0
	nbRank2 = 0
	nbRank3 = 0
	for elt in dico.keys():
		if dico[elt][1] != 'No' or  dico[elt][2] != 'No' or dico[elt][3] != 'No' :
			if dico[elt][0] == 1 :
				nbRank1 += 1
			if dico[elt][0] == 2 :
				nbRank2 += 1
			if dico[elt][0] == 3 :
				nbRank3 += 1
				
######################### Create output file ########################################

	with open('%s%s_compare.txt'%(outDir,Id),'w') as f
		if signalp == 'None':
			f.write('There is no output file for SignalP by default the SignalP column contains only "No"\n')
		if targetp == 'None':
			f.write('There is no output file for TargetP by default the TargetP column contains only "No"\n')
		if phobius == 'None':
			f.write('There is no output file for Phobius by default the Phobius column contains only "No"\n')
		f.write('Rank : Reliability class, from 1 to 3, where 1 indicates the strongest prediction.\n\n')


		f.write(f'Prediction Info :\n\t- Total : {len(listeGene)} gene ({len(dico)} transcrit)\n'
				f'\t- Predict by SignalP : {len(GeneSignalP)} gene ({len(listesignalp)} transcrit)\n'
				f'\t- Predict by targetp {len(GeneTargetP)} gene ({len(listetargetp)} transcrit)\n'
				f'\t- Predict by phobius : {len(GenePhobius)} gene ({len(listephobius)} transcrit)\n\n'
				f'Info prédiction :\n\t- Number of  rank1 (Predict by all tools) : {len(listeRank1)} gene ({nbRank1} transcrit)\n\n'
				f'\t- Number of rank2 (Predict by 2 tools) : {len(listePhTar)+len(listePhSign)+len(listeSignTar)} gene ({nbRank2} transcrit)\n'
				f'\t\t- Predict by signalP and targetP : {len(listeSignTar)} gene ({nbSignTar} transcrit)\n'
				f'\t\t- Predict by signalP and phobius : {len(listePhSign)} gene ({nbPhSign} transcrit)\n'
				f'\t\t- Predict by targetP and phobius : {len(listePhTar)} gene ({nbPhTar} transcrit)\n\n'
				f'\t- Number of rank3 (Predict by only 1 tools) : {len(listePh)+len(listeSign)+len(listeTar)} gene ({nbRank3} transcrit)\n'
				f'\t\t- Predict by signalP only : {len(listeSign)} gene ({nbSign} transcrit)\n'
				f'\t\t- Predict by targetP only : {len(listeTar)} gene ({nbTar} transcrit)\n'
				f'\t\t- Predict by phobius only: {len(listePh)} gene ({nbPh} transcrit)\n\n')


		lengthID = []
		for elt in dico.keys() :
			lengthID.append(len(elt))
		lenID = max(lengthID)

		f.write(' -%-12s---%5s---%7s---%7s---%7s- \n'%("-".center(lenID,'-'),"-"*5,"-"*7,"-"*7,"-"*7))
		f.write('| %s | %5s | %7s | %7s | %7s |\n'%("Gene_id".center(lenID),"Rank".center(5),"SignalP","TargetP","Phobius"))
		f.write('| %-12s---%5s---%7s---%7s---%7s |\n'%("-".center(lenID,'-'),"-"*5,"-"*7,"-"*7,"-"*7))
		for elt  in sorted(dico.keys(), key=sort_human):
			if dico[elt][1] != 'No' or  dico[elt][2] != 'No' or dico[elt][3] != 'No' :
				f.write('| %s | %5s | %7s | %7s | %7s |\n'%(elt.center(lenID),str(dico[elt][0]).center(5) ,dico[elt][1].center(7),dico[elt][2].center(7),dico[elt][3].center(7)))
				listeFinal.append([dico[elt][0],dico[elt][1]])

	dico_fasta = fasta2dict(fasta)
	with open(f'{outDir}{Id}_secreted_1.fasta','w') as f :
		for idSeq in sorted(dico.keys(), key=sort_human):
				if dico[idSeq][0] <= rankMini :
					seqObj = dico_fasta[idSeq].seq
					record = SeqRecord(seqObj,id=idSeq,name=idSeq, description= dico_fasta[idSeq].description + ' | rank_'+str(dico[idSeq][0]))
					SeqIO.write(record,f, "fasta")

	
	
############## summary message #######################

	print(form('\n-----------------------------------------------------------------------------------------------------------------------','red','bold'))
	print(form('Execution summary:\n','green',['bold','underline']))
	print('\n\tInput : \n\t\t- '+ fasta+'\n\t\t-'+signalp+'\n\t\t-'+targetp+'\n\t\t-'+phobius)
	print('\n\tOutput :')
	print('\t\t - Résultat des prédictions des secretomes : '+outDir)
	print('\n\tRank mini : %s\n'%rankMini)
	print('\tInfo prediction :\n\t\t- Total transcrit : %s\n\t\t- Gene predit par signalp %s\n\t\t -Gene predit par targetp %s\n\t\t -Gene predit par phobius %s\n\n\tInfo prediction :\n\t\t- Nombre de transcrit rank1 : %s\n\t\t- Nombre de transcrit rank2 : %s\n\t\t- Nombre de transcrit rank3 : %s'%(len(dico),len(listesignalp),len(listetargetp),len(listephobius),nbRank1,nbRank2,nbRank3))
	print(form('----------------------------------------------------------------------------------------------------------------------','red','bold'))

		
	
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))



