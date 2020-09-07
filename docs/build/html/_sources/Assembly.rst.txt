
Assembly Script
===============

To assemble our pair-end data, we used the ABySS tool with different k parameters (kmer) for each isolate. This assembly was launched with ABYSS_launch.py script.

ABYSS_launch
------------

This script does not directly launch ABySS but creates output directory and bash scripts to launch all assembly in parallel. This script is used on cirad HPC `link <https://bioinfo-agap.cirad.fr/>`_\ , so the environment (module load command) and sge submission works only on the cirad cluster. Please adapt this script to your machine or use the new pipeline which uses the singularity containers and configuration cluster file. 

- \ `ABySS <https://github.com/bcgsc/abyss>`_\ 
- \ `Python >=3.7 <https://www.python.org/downloads/>`_\ 

Arguments take by ABYSS_launch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


- **-d, --directory** (type : string) : path of directory that contains all the fastq files (with "_R1.fastq.gz" and
  "_R2.fastq.gz" extension ) which must be assembled
- **-o, --out** (type : string) : path of output directory for all output file (script, log, stat and assembly)

**Exemple:**

.. code-block::

   ./ABYSS_launch.py -d /homedir/user/work/data/ -o /homedir/user/work/result/

formatFastaName
---------------

The fomatFastaName.py is used to renamed correctly the scaffold of assembly in function of the length ( exemple : scaffold_1 is the longer scaffold et scaffold_2 the second more long) and filter scaffold in function of the length.

Mandatory installation
~~~~~~~~~~~~~~~~~~~~~~


- \ `Python >=3.7 <https://www.python.org/downloads/>`_\ 

Arguments take by formatFastaName
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


- **-f, --fasta** (type : string) : path of the fasta to rename and filter
- **-o, --out** (type : string) : path of the output file
- **-k, --keep** (type : string) : choice keep sequences size greater than -l (g/greater) or keep lower (l/lower
  ) than **-l, --len** parameter
- **-l, --len** (type : integer) : scaffold length cutoff

**Exemple:**

.. code-block::

   # Keep sequences greater than 1000
   format_fasta_name.py -f sequences.fasta -l 1000 -o sequence_Sup1000.fasta -k g

   # Keep sequences lower than 1000
   format_fasta_name.py -f sequences.fasta -l 1000 -o sequence_Inf1000.fasta -k l
