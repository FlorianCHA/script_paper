
Secretome Selection
===================

For select secretome from protein fasta of many sample, we use secretome_Piepline. This script is used on `cirad HPC <https://bioinfo-agap.cirad.fr/>`_\ , so the environment (module load command) and sge submission works only on the cirad cluster. This script does not directly launch ABySS but creates output directory and bash scripts to launch all assembly in parallel.

Secretome_Pipeline
------------------

This program is used to predict secretome. This program uses :


#. SignalP, TagetP, Phobius for detect the presence of signal peptide cleavage sites.
#. ComparaisonSecretome script to retrieve and compare information from the secretome prediction tools.
#. Selection_TMHMM script to select, from TMHMM ouput, only protein with no transmembrane domain or only one
   transmembrane in 40 first aa.
#. EliminateREmotif script to eliminate, from ps_scan ouput, the protein with a RE retention motif.

All step hasn't launch by the script, at the end of the script a message give you a command line to enter in your
 terminal for launch all job

Mandatory installation
~~~~~~~~~~~~~~~~~~~~~~


- \ `SignalP == 4.1 <http://www.cbs.dtu.dk/services/SignalP-4.1/>`_\ 
- \ `targetP <http://www.cbs.dtu.dk/services/TargetP/>`_\ 
- \ `phobius <http://phobius.sbc.su.se/data.html>`_\ 
- \ `PS-scan <https://prosite.expasy.org/scanprosite/>`_\ 
- \ `TMHMM <http://www.cbs.dtu.dk/services/TMHMM/>`_\ 
- \ `Python >=3.7 <https://www.python.org/downloads/>`_\ 

Arguments take by secretome_Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


- **-f, --file** (type : string) : path of fasta files that contains all the protein
- **-o, --outdir** (type : string) :  path of output directory for all output file (script, log, stat and assembly)
- **-p, --prosite** (type : string) : path of prosite.dat file. You can upload the file **ftp://ftp.expasy.org/databases/prosite/prosite.dat**
- **-fo, --force** (type : none) : Force the script to remove output data

**Exemple:**

.. code-block::

   secretome_Pipeline.py -d /homedir/user/work/data/ -o /homedir/user/work/result/ -p /homedir/user/work/prosite.dat

Comparaison_secretome
---------------------

This program is used to retrieve and compare information from the output of secretome prediction tools (signalP, targetP and Phobius) and create a summary table. This program is used by secretome_Pipeline.

Mandatory installation
~~~~~~~~~~~~~~~~~~~~~~


- \ `Python >=3.7 <https://www.python.org/downloads/>`_\

Arguments take by comparaison_secretome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


- **-s, --signalp** (type : string) : path of the signalP output file
- **-t, --targetp** (type : string) : path of the targetp output file
- **-p, --phobius** (type : string) : path of the phobius output file
- **-f, --fasta** (type : string) : path of the input fasta file used in singalP, targetp and phobius
- **-r, --rank** (type : int) :  rank mini for selection by default rank = 3. 

  * rank 1 = protein secreted detect by all tools
  * rank 2 = protein secreted detect by two tools,
  * rank 3 : protein secreted detect by one tool 

- **-o, --outdir** (type : string) :  path of output directory for all output file

**Exemple:**

.. code-block::

   # Comparaison between all tools
   comparaisonSecretome.py  --signalp /homedir/user/work/data/output_signalp.txt --targetp /homedir/user/work/data/output_targetp.txt --phobius /homedir/user/work/data/output_phobius.txt -f  /homedir/user/work/protein.fasta -o /homedir/user/work/result/

   # Comparaison and selection with two tools
   comparaisonSecretome.py  --signalp /homedir/user/work/data/output_signalp.txt --targetp /homedir/user/work/data/output_targetp.txt -f  /homedir/user/work/protein.fasta -o /homedir/user/work/result/

   #Comparaison and selection with only one tools
   comparaisonSecretome.py  --signalp /homedir/user/work/data/output_signalp.txt -f  /homedir/user/work/protein.fasta -o /homedir/user/work/result/

EliminateREmotif
----------------

This program is used to eliminate sequence from the result of PS-scan with PS00014 pattern (RE target Motif). This program is used by secretome_Pipeline.

Mandatory installation
~~~~~~~~~~~~~~~~~~~~~~


- \ `Python >=3.7 <https://www.python.org/downloads/>`_\ 

Arguments take by eliminateREmotif
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


- **-p, --ps_scan** (type : string) : path of the ps_scan output file
- **-f, --fasta** (type : string) : path of the fasta file which ps_scan has been proceed
- **-o, --outdir** (type : string) :  path of output directory for all output file

**Exemple:**

.. code-block::

   eliminateREmotif.py  -p /homedir/user/work/data/TMHMM_result.txt -o /homedir/user/work/result/ -f  /homedir/user/work/result/protein_filter.fasta

SelectionTMHMM
--------------

This program is used to selected the proteins having no transmembrane domain or only one in the 40 first amino acid. This program is used by secretome_Pipeline.

Arguments take by selectionTMHMM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


- **-t, --TMHMM** (type : string) : path of the TMHMM output file
- **-f, --fasta** (type : string) : path of the fasta file which TMHMM has been proceed
- **-o, --outdir** (type : string) :  path of output directory for all output file

**Exemple:**

.. code-block::

   selection_TMHMM.py  -t /homedir/user/work/data/TMHMM_result.txt -o /homedir/user/work/result/ -f  /homedir/user/work/result/protein_filter.fasta
