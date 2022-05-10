# gff_to_multigene_fasta
A script to make a multigene fasta file from a mitogenome fasta file and gff file from the GenBank, created for personal use.
If we have a .gff3 file and a fasta file in a folder, we can make a multigene fasta file from these two files.

Often it is essential to get separate gene alignments for phylogenetics, and it is quite possible that we are using a mitogenomic dataset. Some species' data will be available only as a complete mitogenome file in the GenBank. This script reads a complete mitogenome file and make a multigene fasta file for you.

Please open the file to see the code.

provide input when prompted.

_PS: It creates a lot of junk files in the folder where you work from, please delete them since they are of no use probably._

**Usage**

>python /path/to/where/the/script/is/stored/gtt_to_multigene_fasta.py

input your fasta file name: eg. **trial.fasta**

Resultant file with multigenes will be - trial_multigene.fasta (or ***_multigene.fasta**).

**Absolute Requirement**

The folder where you work from should have a file with extension _.gff3_ for it to work.

