#!/usr/bin/python
import csv
from operator import itemgetter
import re
from Bio import SeqIO
import os

file_in = raw_input('\n\n please give the name of the NCBI mitogenome fasta file: ')
for file in os.listdir("."):
	if file.endswith("gff3"):
		gff_f = file
gff_fl=gff_f
adc = file_in.split('.')
abc = file_in.split('.')
mg="_multigene"
mg_n = adc[0]+mg
sgn = "_single_line.fasta"
out_file = adc[0]+sgn
ofile=open(out_file,'w')

#for reverse complementing - strand sequences from https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python/25189373
alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases

#three new files are written in the current folder - getting gene coordinates from gff file and preparing a table of it in a file, multiline to single line fasta and in a new file write only the sequence without header
def get_gene_cords(gff_file, fasta_file):
	gene_l = []
	gene_r = []
	gene_n = []
	gene_names = []
	polarity = []
	with open(gff_fl) as g:
    		rd = csv.reader(g, delimiter="\t", quotechar='"')
    		for row in rd:
			if len(row) < 1:
				continue
			elif len(row) > 1:
				if row[2] == 'gene':
					gene_l.append(row[3])
					gene_r.append(row[4])
					polarity.append(row[6])
					gene_n.append(row[8])
    		for item in gene_n:
    			g_n = re.findall("gene\-(.*?)\;", item)
			gene_names.append(g_n)
	rows = zip(gene_l,gene_r,gene_names,polarity)

	with open("cds_table.csv", "w") as f: #write the coordinates and gene names into a csv file
		writer = csv.writer(f)
		for row in rows:
			writer.writerow(row)

	with open(file_in, 'r') as myfile: #multi line fasta to single line fasta
		for record in SeqIO.parse(myfile, "fasta"):
			sequence = str(record.seq)
			ofile.write('>'+record.id+'\n'+sequence+'\n')

	with open(out_file) as f:
		next(f)		#get only the second line of the fasta file
		for line in f:
			seq_line=line
			seq_only_file=open('seq_only.txt','w')
			seq_only_file.write(line) #write the sequence only in a new file to read in the next step
			seq_only_file.close()
	
	with open('cds_table.csv', 'rb') as f:
    		reader = csv.reader(f)
    		for row in reader:
			left=[]
			left.append(row[0])
			right=[]
			right.append(row[1])
			name=[]
			name.append(row[2])
			b = 0
			while b < len(left):	
				lc = left[b]
				rc = right[b]
				left_cord = int(lc)-1
				right_cord = int(rc)-1
				fasta=">"
				header=name[b]
				fasta_line=fasta+header
				with open('seq_only.txt','r') as seq_line:
					sequence=seq_line.readline()
					#print fasta_line+'\n'+sequence[left_cord:right_cord]+'\n'
					mgn="_multigene.fasta"
					mg_nn=abc[0]+mgn
					genes=open(mg_nn,'a')
					if row[3] == "-":
						rcp = int(rc)
						minus_seq=sequence[left_cord:rcp]
						m_seq = reverse_complement(minus_seq)
						genes.write(fasta_line+'\n'+m_seq+'\n')
					elif row[3] == "+":
						genes.write(fasta_line+'\n'+sequence[left_cord:right_cord]+'\n')
						genes.close()
				b = b+1
	for file in os.listdir("."): #prepare a new final fasta file without the square brackets as above
		if file.endswith("_multigene.fasta"):
			#print file
			fl = file
			#print fl
			with open(fl,'r') as ch_fasta:
				chf=ch_fasta.readlines()
				for line in chf:
					changed_lines=re.sub(r"\['|\']", r"", line)
					#print changed_lines
					mgnf="_final_multigene.fasta"
					mg_nnf=abc[0]+mgnf
					final=open(mg_nnf,'a')
					final.write(changed_lines)
					final.close()
get_gene_cords(gff_fl, file_in)
