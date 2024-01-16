#!/usr/bin/env python
import sys
import TraDIS_simu_funct
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import random

#Define parameters
Mariner = 'GCCAACCTGT'
insert_counter = 0
insert_lib ={}
insert_positions = []
selected_insert_positions = []
fasta_lib={}

#import and load genomic data file & insertion no into python & check
if len(sys.argv)<3:
	print('invalid command line argument: include a genbank file and a insertion number')
	sys.exit(0)
genbank_file=sys.argv[1] 
insertions=int(sys.argv[2])

#parse genomic data and essential genes from the file
seq_id,genome_data=TraDIS_simu_funct.file_parse_genome(genbank_file)
essential_genes=TraDIS_simu_funct.essential_gene_parse(genbank_file)

#extend genome 
three_prime_extension=genome_data[0:40]
five_prime_extension=genome_data[-40:]
ext_genome=five_prime_extension+genome_data+three_prime_extension

#insert transposon at TA sites
rand_insert_positions=TraDIS_simu_funct.transposon_insert(genome_data,insert_positions,selected_insert_positions,insertions)

#remove insertions that have inserted into essential genes
rand_insert_positions=TraDIS_simu_funct.essential_gene_insert_check(essential_genes,rand_insert_positions)

#randomly select transposons up to specified insertion number
insert_no=0
insert_prob=random.randrange(10) 
for insert in rand_insert_positions:
	insert_prob=random.randrange(10)
	insert_no+=1
	if insert_prob < 5:
		insert_lib=TraDIS_simu_funct.transposon_seq_down(insert,ext_genome,insert_lib,insert_no)
	if insert_prob >=5:
		insert_lib=TraDIS_simu_funct.transposon_seq_up(insert,ext_genome,insert_lib,insert_no)

#write sequence output in fasta format
for key,seq in insert_lib.items():
	fasta_lib=TraDIS_simu_funct.fasta_format(fasta_lib,seq,key,Mariner,seq_id)
	
#write fasta formatted outputs to file in fasta format
TraDIS_simu_funct.fasta_file_write(fasta_lib,seq_id,insertions)

#write essential gene locus tags to file
TraDIS_simu_funct.essential_gene_out(essential_genes,seq_id,insertions) 
