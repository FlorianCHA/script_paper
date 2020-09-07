#!/usr/local/bioinfo/python/3.6.4/bin/python
# -*- coding: utf-8 -*-
# @package module_Flo.py
# @author Florian Charriat
# __docformat__ = "restructuredtext en"

"""
    The module_Flo module
    =====================

    :author: CHARRIAT Florian\n
    :contact: florian.charriat@inra.fr\n
    :date: 21/03/2018\n
    :version: 0.1\n

    Use it to import very handy functions.

    Example:

    >>> from module_Flo import createDir
    >>> createDir('resultat')

"""
##################################################
## Modules
##################################################
## Python modules
import argparse, os, glob, re, itertools, numpy as np,  collections, statistics

from math import *
#import statistics
## BioPython
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
## for parse
from collections import namedtuple
import gzip
import urllib


# import statistics


####### FUNCTION ################

def is_number(s):
    """
    This function test if "s" is a float
    :Paramerters:
        s (float) : the float to test
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


################################ directory function ##################################################"

def createDir(Listedirectory):
    '''
    Allows you to check if a folder exists, if not, the folder will be created.

    :Parameters:
         Listedirectory (list) :  list of directory to create
    '''

    if type(Listedirectory) != list:
        if not os.path.exists(Listedirectory):
            os.makedirs(Listedirectory)

    else:
        for directory in Listedirectory:
            if not os.path.exists(directory):
                os.makedirs(directory)
    return

def verifDir(directory, check=False):
    '''
    Allows to format the folder path to be used in a script, the function checks if there is a '/'
    at the end of the path, otherwise it adds it. The function can also verify that a directory exists.

    :Parameters:
         directory (str) :   Path du dossier
         check (bool) : if check = True, the function verify that the directory exists

    '''
    if directory.endswith('/') == False:
        directory = directory + '/'
    if check:
        if os.path.isdir(directory):
            return directory
        else:
            raise ValueError(
                form("ERROR the directory '%s' is not valid path, please check if your directory exists" % directory,
                     'red', 'bold'))
    else:
        return directory


################################## Fonction fichier ############################################"

def verifFichier(file):
    '''
    Allows you to check if a file exists.
    :Parameters:
        file (str) : Path du fichier
    '''
    if os.path.exists(file):
        return
    else:
        raise ValueError(
            form("ERROR the file '%s' doesn't exist, please check if your files exists" % file, 'red', 'bold'))


##################################### Fonction fichier fasta/fastq #################################################"

def isFasta(file):
    '''Used to check if a file is in fasta format, return True if the file is in fasta format
    :Parameters:
        file (str) :  Path du fichier
    '''
    if file.endswith('.fasta') or file.endswith('.fa') or file.endswith('.fasta.gz') or file.endswith(
            '.fa.gz') or file.endswith('.fna'):
        return True
    else:
        return False

def isFastq(file):
    '''Used to check if a file is in fastq format, return True if the file is in fasta format
    :Parameters:
        file (str) :  Path of file
    '''
    if file.endswith('.fastq') or file.endswith('.fq') or file.endswith('.fastq.gz') or file.endswith(
            '.fq.gz'):
        return True
    else:
        return False

def recupId(file):
    '''Allows you to retrieve the name of the file without the fasta or fastq extension or what is after the '_'=
    :Parameters:
         file (str) : Path of file
    '''
    # Traitement pour fichier fasta
    file = file.replace('.fasta.gz', '')
    file = file.replace('.fa.gz', '')
    file = file.replace('.fasta', '')
    file = file.replace('.fa', '')

    # Traitement pour fichier fastq
    file = file.replace('.fastq.gz', '')
    file = file.replace('.fq.gz', '')
    file = file.replace('.fastq', '')
    file = file.replace('.fq', '')

    # Si id avec '_' garde seulement le premier identifiant
    file = file.split('_')[0]
    return file

def fasta2dict(file):
    """
    Function that take a file name (fasta), and return a dictionnary of sequence
    :Parameters:
        file (str) :  Path of file
    """
    with open(file, "r") as fastaFile:
        return SeqIO.to_dict(SeqIO.parse(fastaFile, "fasta"))

def lenSeq2dict(file):
    """
    Function that take a file name (fasta), and return a dictionnary with length of sequence
    :Parameters:
        file (str) : Path of fasta file
    """
    dico_length = {}
    dico_fasta = fasta2dict(file)
    for gene in sorted(dico_fasta.keys(), key=sort_human):
        if dico_fasta[gene].id not in dico_length:
            lenseq = len(dico_fasta[gene].seq)
            dico_fasta[gene]=int(lenseq)
    return dico_fasta

def header_script(txt, version):
    """
    Function that create a header for script file
    :Parameters:
        txt (str) : name to the script
        version (str/int) : version of the script
    """
    txt_format = f'|     Welcome in {txt} (Version {version})     |'
    print(form('\n\t' + '-' * len(txt_format), 'yellow', 'bold'))
    print(form('\t' + txt_format, 'yellow', 'bold'))
    print(form('\t' + '-' * len(txt_format) + '\n', 'yellow', 'bold'))

def footer_script() :
    """
    Function that create a footer for script file
    """
    txt_format = f'|          End of execution          |'
    print(form('\n\t' + '-' * len(txt_format), 'yellow', 'bold'))
    print(form('\t' + txt_format, 'yellow', 'bold'))
    print(form('\t' + '-' * len(txt_format) + '\n', 'yellow', 'bold'))

#################################### Fontion formatage texte ################################################

def sort_human(s, _nsre=re.compile('([0-9]+)')):
    """
    Sort the list in the way that humans expect, use list.sort(key=sort_human) or sorted(list, key=sort_human)).
    """
    try:
        return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]
    except TypeError:
        if not isinstance(s, int):
            print("WARNNING MODULES_SEB::sort_human : List %s value not understand so don't sort \n" % s)

    return s


def form(text, col='white', type='none'):
    '''
    Used to format the texts displayed on the terminal.
    :Parameters:
         text (str) : Text to format
         col (str) : The desired color between the colors red, green, yellow, orange, blue and purple
         type (str/list)  : bold, underline, blind et highligth
    '''
    W = '\033[0'  # white (normal)
    R = '\033[31'  # red
    G = '\033[32'  # green
    Y = '\033[33'  # yellow
    O = '\033[33'  # orange
    B = '\033[34'  # blue
    P = '\033[35'  # purple
    end = '\033[0m'  # white (normal)
    Bold = ';1'
    underline = ';4'
    blind = ';5'
    highlight = ';7'
    text = 'm' + text
    if 'bold' in type:
        text = Bold + text
    if 'underline' in type:
        text = underline + text
    if 'highlight' in type:
        text = blind + text
    if 'highlight' in type:
        text = highlight + text
    if col == 'red':
        return R + text + end
    elif col == 'white':
        return W + text + end
    elif col == 'green':
        return G + text + end
    elif col == 'yellow':
        return Y + text + end
    elif col == 'orange':
        return O + text + end
    elif col == 'blue':
        return B + text + end
    elif col == 'purple':
        return P + text + end

def parseBlast6(path):
    '''
    This function take as input a blast result file parth at format 6 (outfmt 6 option in blast command).
    In output, the function give a dictionnay which contains all hits sorted by sequence Subject.
    :Parameters:
         path (str) :  Path of blast result file
    '''
    dico = {}
    with open(path, 'r') as blast_file:
        for line in blast_file:
            query, subject, identity, len_aln, _, _, start_query, end_query, stfzart_subject, end_subject, evalue, bitscore = line.split()
            if subject not in dico.keys():
                dico[subject] = {"identity": [float(identity)], "len_query": [int(len_aln)],
                                 "len_subject": [int(abs(int(end_subject) - int(start_subject)))],
                                 "start_query": [int(start_query)], "end_query": [int(end_query)],
                                 "start_subject": [int(start_subject)], "end_subject": [int(end_subject)],
                                 "evalue": [float(evalue)], "bitscore": [float(bitscore)], 'type': 'Single-Hit'}
            else:
                dico[subject]["identity"].append(float(identity))
                dico[subject]["len_query"].append(int(len_aln))
                dico[subject]["len_subject"].append(int(abs(int(end_subject) - int(start_subject))))
                dico[subject]["start_query"].append(int(start_query))
                dico[subject]["end_query"].append(int(end_query))
                dico[subject]["start_subject"].append(int(start_subject))
                dico[subject]["end_query"].append(int(end_query))
                dico[subject]["start_subject"].append(int(start_subject))
                dico[subject]["end_subject"].append(int(end_subject))
                dico[subject]["evalue"].append(float(evalue))
                dico[subject]["bitscore"].append(float(bitscore))
                dico[subject]["type"] = 'Multi-Hit'

    return (dico)

def parseBlast(path):
    '''
    This function take as input a blast result file parth at normal format.
    In output, the function give a dictionnay which contains all hits sorted by sequence Subject.
    :Parameters:
         path (str) : Path of blast result file
    '''
    dico = {}
    start_aln = True
    subject ='None'
    start_query = True
    with open(path, 'r') as blast_file:
        for line in blast_file:

            if line[0:7] == 'Length=' and start_query == True :
                len_aln = line.replace('\n','').replace('Length=','').strip()
                start_query  = False
            if line[0] == ">":
                end_aln = False
                if subject != 'None' :
                    if subject not in dico.keys() :
                        dico[subject] = {"identity": [float(identity)], "len_query": [int(len_aln)],
                                         "len_subject": [int(abs(int(end_subject) - int(start_subject)))],
                                         "start_query": [int(start_query)], "end_query": [int(end_query)],
                                         "start_subject": [int(start_subject)], "end_subject": [int(end_subject)],
                                         "evalue": [float(evalue)], "bitscore": [int(bitscore)], 'type': 'Single-Hit',
                                         "codon_start": codon_start, 'codon_stop': codon_stop,'strand': [strand]}
                    else:
                        dico[subject]["identity"].append(float(identity))
                        dico[subject]["strand"].append(strand)
                        dico[subject]["len_subject"].append(int(abs(int(end_subject) - int(start_subject))))
                        dico[subject]["start_query"].append(int(start_query))
                        dico[subject]["end_query"].append(int(end_query))
                        dico[subject]["start_subject"].append(int(start_subject))
                        dico[subject]["end_subject"].append(int(end_subject))
                        dico[subject]["evalue"].append(float(evalue))
                        dico[subject]["bitscore"].append(int(bitscore))
                        dico[subject]["type"] = 'Multi-Hit'
                        if strand == 'Minus':
                            dico[subject]["codon_start"] = codon_start
                        else:
                            dico[subject]["codon_stop"] = codon_stop
                subject = line.replace('>', '').replace('\n', '').split()[0].strip()
            elif 'Score =' in line and 'Expect =' in line and end_aln == False :
                evalue = line.split('Expect =')[-1].replace('\n','').strip()
                bitscore = float(line.split('bits')[0].replace('Score =', '').strip())
                start_aln = True
            elif 'Identities =' in line:
                identity = line.split('(')[1].split(')')[0].replace('%','')
            elif 'Query ' in line and start_aln == True:
                start_query = line.split()[1].strip()
                end_query = line.split()[3].replace('\n', '').strip()
            elif 'Query ' in line and start_aln == False:
                end_query = line.split()[3].replace('\n', '').strip()
            elif 'Sbjct' in line and start_aln == True:
                start_subject = line.split()[1].strip()
                end_subject = line.split()[3].replace('\n', '').strip()
                codon_start = line.split()[2].strip()[0:3]
                codon_stop = line.split()[2].strip()[-3:]
                start_aln = False
                end_aln = True
            elif 'Sbjct' in line and start_aln == False:
                end_subject = line.split()[3].replace('\n', '').strip()
                codon_stop = line.split()[2].strip()[-3:]
                end_aln = True
            elif 'Strand=' in line:
                strand = line.split('/')[-1].replace('\n', '').strip()
            elif 'Score =' in line and 'Expect =' in line and end_aln == True:
                evalue = line.split('Expect =')[-1].replace('\n','').strip()
                bitscore = float(line.split('bits')[0].replace('Score =', '').strip())
                start_aln = True
                end_aln = False
                if subject not in dico.keys() :
                    dico[subject] = {"identity": [float(identity)], "len_query": [int(len_aln)],
                                     "len_subject": [int(abs(int(end_subject) - int(start_subject)))],
                                     "start_query": [int(start_query)], "end_query": [int(end_query)],
                                     "start_subject": [int(start_subject)], "end_subject": [int(end_subject)],
                                     "evalue": [float(evalue)], "bitscore": [int(bitscore)], 'type': 'Single-Hit',
                                     "codon_start": codon_start, 'codon_stop': codon_stop,'strand': [strand]}
                else:

                    dico[subject]["identity"].append(float(identity))
                    dico[subject]["len_subject"].append(int(abs(int(end_subject) - int(start_subject))))
                    dico[subject]["start_query"].append(int(start_query))
                    dico[subject]["end_query"].append(int(end_query))
                    dico[subject]["start_subject"].append(int(start_subject))
                    dico[subject]["end_subject"].append(int(end_subject))
                    dico[subject]["evalue"].append(float(evalue))
                    dico[subject]["bitscore"].append(int(bitscore))
                    dico[subject]["strand"].append(strand)
                    dico[subject]["type"] = 'Multi-Hit'
                    if strand == 'Minus':
                        dico[subject]["codon_start"] = codon_start
                    else:
                        dico[subject]["codon_stop"] = codon_stop

    return (dico)

def selectBlast(path, indentity_min, max_intron):
    """
    This function use the parseBlast function for parse the blast result file give with the path option.
    In output, tje function give a list of sequence which is select with the parameters : the minimum of identity and the max length of intron.
    :Parameters:
         path (str) : Path of blast result file
         indentity_min (int) : Minimum of identity to select sequence
         max_intron (int) : Maximum length of intron to select sequence
    """
    dico_blast = parseBlast(path)
    select = []
    for elt in dico_blast.keys():
        hit = dico_blast[elt]
        if hit['type'] == 'Single-Hit':
            if hit["identity"][0] >= indentity_min and hit['len_query'][0] * 0.90 < hit['len_subject'][0] < \
                    hit['len_query'][0] * 1.1:
                select.append(elt)
        if hit['type'] == 'Multi-Hit':
            nbHit = len(hit['identity'])
            for i in range(0, (nbHit - 1)):
                len_intron = hit['end_subject'][i] - hit['start_subject'][i + 1]
                if len_intron <= max_intron:
                    intron = True
                else:
                    intron = False
                    break
            if intron and statistics.mean(hit["identity"]) >= indentity_min and sum(hit['len_query']) * 0.90 < sum(
                    hit['len_subject']) < sum(hit['len_query']) * 1.1:
                select.append(elt)
    return (select)

def functionSens(pos1, pos2):
    """
    Look at the positions of a genomic element to give the meaning of this element,
    and put the positions in the correct order (ascending order)

    :Parameters:
         pos1 (int) : First position of element
         pos2 (int) : Second position of element
    """
    if pos1 > pos2:
        sens = '-'
        start = pos2
        end = pos1
    if pos1 < pos2:
        sens = '+'
        start = pos1
        end = pos2
    return start, end, sens

def indexEgale(liste, target):
    '''
    Function that returns the position of all elements equal to the target in a given list
    :Parameters:
         liste (list) :  The list where you search the target
         target (str)  Target to serach in the list
    '''
    index = []
    for i, e in enumerate(liste):
        if e == target:
            index.append(i)
    return index


def indexDif(liste, target):
    '''
    Function which allows to return the position of all the elements different from the target str in a given list
    :Parameters:
         liste (list) :  The list where you search the target
         target (str)  Target to remove of the list
    '''
    index = []
    for i, e in enumerate(liste):
        if e != target:
            index.append(i)
    return index

def comparaisonListe(list1, list2):
    """
    Allows you to retrieve the elements in common between the two lists (duplicates are eliminated)
    :Parameters:
         list1 (list) :  The first list to compare
         list2 (list) : The second list to compare
    """
    list = []
    list1 = set(list1)
    list2 = set(list2)
    for elt in list1:
        if elt in list2:
            list.append(elt)
    return (list)

def IdisIn(list, text):
    """
    Allows you to check if an element of the list is in the text
    :Parameters:
         list (list) :  The list which contain all elt to test
         text (list) : The target text
    """
    list = set(list)
    for elt in list:
        if elt in text:
            return True, elt
    return False, None

class genpop :
    """
    Object making it possible to perform certain population genetic analyzes (For the moment, only nucleotide diversity is used    """
    def __init__(self,aln):
        self.aln = aln
        self.Pi = 'Not yet calculate, please use the ".calculate_statistics" function for calculate'
        self.Tajima_D = 'Not yet calculate, please use the ".calculate_statistics" function for calculate'
        self.lseff = 'Not yet calculate, please use the ".calculate_statistics" function for calculate'
        self.nseff ='Not yet calculate, please use the ".calculate_statistics" function for calculate'
        self.nb_haplotype = 'Not yet calculate, please use the ".calculate_statistics" function for calculate'
        self.haplotype = 'Not yet calculate, please use the ".calculate_statistics" function for calculate'

    def make_aln(self, MISS):
        """
        Create an Align object
        """
        liste_dna = list()
        for id in self.aln .keys():
            record_dna = SeqRecord(Seq(str(self.aln [id].seq).replace('-', 'N').upper()), id=id, name=id,
                                   description=self.aln [id].description)
            liste_dna.append(record_dna)
        aln_nucl= MultipleSeqAlignment(liste_dna)
        nbSeq = len(aln_nucl)
        length_Seq = len(aln_nucl[1, :].seq)
        for i in range(0, length_Seq):
            str_columns = aln_nucl[:, i].upper()
            if ((str_columns.count('-')+str_columns.count('N')) / nbSeq) > MISS:
                continue
            else :
                try:
                    edited += aln_nucl[:, i:i+1]
                except:
                    edited = aln_nucl[:, i:i+1]
        try:
            edited
        except NameError:
            edited = ''
        self.aln = edited


    def make_alnS(self, MISS):
        """
        Create an Align with only Synonymous codon sites
        """
        liste_dna = list()
        liste_prot = list()
        for id in self.aln .keys():
            record_dna = SeqRecord(Seq(str(self.aln [id].seq).replace('-', 'N').upper()), id=id, name=id,
                                   description=self.aln [id].description)
            liste_dna.append(record_dna)
            record = SeqRecord(Seq(str(self.aln [id].seq).replace('-', 'N')).translate(), id=id, name=id,
                               description=self.aln [id].description)
            liste_prot.append(record)
        aln_prot = MultipleSeqAlignment(liste_prot)
        aln_dna = MultipleSeqAlignment(liste_dna)
        nbSeq = len(aln_prot)
        length_Seq = len(aln_prot[1, :].seq)
        for i in range(0, length_Seq):
            str_columns = aln_prot[:, i].upper()
            if (str_columns.count('X') / nbSeq) > MISS:
                continue
            elif nbSeq == (str_columns.count(str_columns[0]) + str_columns.count('X')):
                if nbSeq != (aln_dna[:, (i * 3)].count(aln_dna[:, (i * 3)][0]) + aln_dna[:, (i * 3)].count('N') + aln_dna[:,(i * 3)].count('-')) or \
                        nbSeq != (aln_dna[:, (i * 3 + 1)].count(aln_dna[:, (i * 3 + 1)][0]) + aln_dna[:, (i * 3 + 1)].count('N') + aln_dna[:, (i * 3 + 1)].count('-')) or \
                        nbSeq != (aln_dna[:, (i * 3 + 2)].count(aln_dna[:, (i * 3 + 2)][0]) + aln_dna[:, (i * 3 + 2)].count('N') + aln_dna[:, (i * 3 + 2)].count('-')):
                    try:
                        edited += aln_dna[:, ((i) * 3):((i + 1) * 3)]
                    except:
                        edited = aln_dna[:, ((i) * 3):((i + 1) * 3)]
        try:
            edited
        except NameError:
            edited = ''
        self.aln = edited

    def make_alnNS(self, MISS):
        """
        Create an Align with only No Synonymous codon sites
        """
        liste_dna = list()
        liste_prot = list()
        for id in self.aln.keys():
            liste_dna.append(self.aln[id])
            record = SeqRecord(Seq(str(self.aln[id].seq).replace('-', 'N')).translate(), id=id, name=id,
                               description=self.aln[id].description)
            liste_prot.append(record)
        aln_prot = MultipleSeqAlignment(liste_prot)
        aln_dna = MultipleSeqAlignment(liste_dna)
        nbSeq = len(aln_prot)
        length_Seq = len(aln_prot[1, :].seq)
        for i in range(0, length_Seq):
            str_columns = aln_prot[:, i].upper()
            if (str_columns.count('X') / nbSeq) > MISS:
                continue
            elif nbSeq != (str_columns.count(str_columns[0]) + str_columns.count('X')):
                try:
                    edited += aln_dna[:, ((i) * 3):((i + 1) * 3)]
                except:
                    edited = aln_dna[:, ((i) * 3):((i + 1) * 3)]
            else :
                    edited = ''

        try:
            edited
        except NameError:
            edited = ''
        self.aln = edited

    def calculate_statistics(self):
        """
        Caculate Pi from make_aln_object
        """
        ########## Calculate diversity nucléotidique (Pi) ########
        if self.aln != '':
            liste = list()
            liste_Pi = list()
            num_polymorphism = list()
            for record in self.aln:
                liste.append(str(record.seq).replace('-','N').upper())
            nb_pair_analysis = 0
            for pair in itertools.combinations(liste, 2):
                nb_pair_analysis += 1
                Pi_pair, numP = numdiffs(pair[0], pair[1], 'N')
                liste_Pi.append(Pi_pair)
                num_polymorphism.append(numP)
            liste_Pi = [x for x in liste_Pi if x != 'NA']
            num_polymorphism = [x for x in num_polymorphism if x != 'NA']
            self.Pi = statistics.mean(liste_Pi)

            ####### Calculate Tajima'D score (D) #############
            nbSeq = len(self.aln)
            length_Seq = len(self.aln[1, :].seq)
            # This part calculate the number of segregating site (S)
            S = 0
            for i in range(0, length_Seq):
                str_columns = self.aln[:, i].upper()
                if nbSeq != (str_columns.count(str_columns[0]) + str_columns.count('N') + str_columns.count('-')) :
                    S += 1
            mean_polymorphisme = sum(num_polymorphism)/nb_pair_analysis
            a1 = sum(1/ d for d in range(1,nbSeq))
            a2 = sum(1/ (d*d) for d in range(1,nbSeq))
            b1 = (nbSeq+1)/(3*(nbSeq-1))
            b2 = 2*(nbSeq*nbSeq +nbSeq+3)/(9*nbSeq*(nbSeq-1))
            c1 = b1 -1/a1
            c2 = b2 - ((nbSeq+2)/(a1*nbSeq)) + a2/(a1*a1)
            e1 = c1/a1
            e2 = c2/(a1*a1+a2)
            D = (mean_polymorphisme - S/a1)/(sqrt(e1*S+e2*S*(S-1)))
            self.Tajima_D = D
            ####### Some statistics #############
            self.lseff = length_Seq
            self.nseff = nbSeq
            self.nb_haplotype = len(set(liste))
            self.haplotype = list(set(liste))
            ### nb haplotype
        else :
            self.Tajima_D = 'None'
            self.lseff = 'None'
            self.nseff = 'None'
            self.nb_haplotype = 'None'
            self.haplotype = 'None'
            self.Pi = 'None'
