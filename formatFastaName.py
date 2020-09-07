#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# @package format_fasta_name.py
# @author Sebastien Ravel

"""
    The format_fasta_name script
    =============================
    :author: Charriat Florian
    :contact: florian.charriat@inra.fr
    :date: 9/03/2018
    :version: 0.1

    Script description
    ------------------

    This Programme filter sequences by length

    Example
    -------

    >>> # Keep sequences greater than 1000
    >>> format_fasta_name.py -f sequences.fasta -l 1000 -o sequence_Sup1000.fasta -k g
    >>> # Keep sequences lower than 1000
    >>> format_fasta_name.py -f sequences.fasta -l 1000 -o sequence_Inf1000.fasta -k l

    Help Programm
    -------------

    optional arguments:
        - \-h, --help
                        show this help message and exit
        - \-v, --version
                        display extractSeqFastaFromLen.py version number and exit

    Input mandatory infos for running:
        - \-f <filename>, --fasta <filename>
                        fasta files
        - \-l <int>, --len <int>
                        lensize cutoff
        - \-o <filename>, --out <filename>
                        name of output file
        - \-k <g/greater/l/lower>, --keep <g/greater/l/lower>
                        choice keep sequences size greater than -l (g/greater) or keep lower (l/lower)

"""


##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
from module_Flo import fasta2dict, lenSeq2dict, header_script,footer_script

## Python modules
import argparse
from time import localtime, strftime
## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Variables Globales
version="0.1"

##################################################
## Main code
##################################################
if __name__ == "__main__":

    # Initializations
    # Parameters recovery
    parser = argparse.ArgumentParser(prog='extractSeqFastaFromLen.py', description='''This Programme filter sequences by length and rename ID''')
    parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
                        'Display extractSeqFastaFromLen.py version number and exit')


    filesreq = parser.add_argument_group('Input mandatory infos for running')
    filesreq.add_argument('-f', '--fasta', metavar="<filename>",type=str, required=True, dest = 'fasta_file', help = 'Fasta files')
    filesreq.add_argument('-o', '--out', metavar="<filename>", required=True, dest = 'output_file', help = 'Name of output')


    files = parser.add_argument_group('Input infos for running with default values')
    files.add_argument('-k', '--keep', metavar="<g/greater/l/lower>",type=str, required=False, dest = 'keep', default ='None', help = 'Choice keep sequences size greater than -l (g/greater) or keep lower (l/lower)')
    files.add_argument('-l', '--len', metavar="<int>",type=int, required=False, dest = 'len_size',default = 500, help = 'Lensize cutoff')
    # Check parameters
    args = parser.parse_args()

    #Welcome message
    print(form("\n\t---------------------------------------------------------",'yellow','bold'))
    print("\t"+form("|",'yellow','bold')+form("       Welcome in mergeBraker_augustus   (Version " + version + ")      ",type='bold')+form("|",'yellow','bold'))
    print(form("\n\t---------------------------------------------------------",'yellow','bold'))

    header_script("format_fasta_name",version)



    # Retrieve all parameters give in input
    fasta_file = os.path.abspath(args.fasta_file)
    output_file = os.path.abspath(args.output_file)
    len_size = args.len_size
    keep = args.keep

     # Open output file for write the new fasta
     with open(output_file, "w") as output_handle :
		# Create a dictionary  for sort the sequence from length
        dicoSize = lenSeq2dict(fasta_file)
		# Create a dico for obtain the sequence from ID
        dicoFasta = fasta2dict(fasta_file)
        filename = fasta_file.split('/')[-1].replace('.fasta','')
        nbTotal = len(dicoFasta.keys())
        count=1
        for ID in sorted(dicoSize.keys(), key=dicoSize.get, reverse=True):
            lenSeq = dicoSize[ID]
            sequence = dicoFasta[ID]
            strain = output_file.split('/')[-1].split('_')[0].replace('.fasta','')
            seqName = f'Scaffold_{count}'
            descrip = f'length={lenSeq}'
			# Check if the sequence have more than 20 pb which not a unknow base ('N')
            if str(sequence.seq).count('N') < (len(str(sequence.seq)) - 20) :
				# Eliminate the sequences which do not correspond to the criterion given in input
                if keep == 'g' and lenSeq >= lenSize or (keep == 'l' and lenSeq <= lenSize):
                        record = SeqRecord(sequence.seq,id=seqName,name=seqName, description=descrip)
                        SeqIO.write(record,output_handle, "fasta")
                        count += 1
                elif keep == 'None' :
                        record = SeqRecord(sequence.seq,id=seqName,name=seqName, description=descrip)
                        SeqIO.write(record,output_handle, "fasta")
                        count += 1
            else :
                print(f'{ID} contains only N, so this scaffolds is removed')


    print('\n\nExecution summary:')

    print(f'-Outputting \nIl y a {count} Sequences\nles sequences sont ajouter dans le fichier {output_file}')

    footer_script()




