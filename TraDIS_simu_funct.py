import re
import random
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

#genmoic data parse function 
def file_parse_genome(gb_file):
	for seq_record in SeqIO.parse(open(gb_file),'genbank'):
		seq_id = seq_record.id
		genome_data = str(seq_record.seq)
	return(seq_id,genome_data)
	
#essentail gene parse function
def essential_gene_parse(gb_file):
	gene_dict = {} 
	rand_genes_dict = {}   
	for rec in SeqIO.parse(open(gb_file), "genbank"):
		if rec.features:
			for feature in rec.features:
				if feature.type == "CDS":
					selected_gene=feature.location
					gene_tag=(feature.qualifiers['locus_tag'])
					locus_id=gene_tag[0]
					gene_dict[locus_id]=selected_gene			
	rand_locus_tags=random.sample(gene_dict.keys(),200)
	for i in rand_locus_tags:
		if i in gene_dict.keys():
			key=i
			value=gene_dict[i]
			rand_genes_dict[key]=value		
	return(rand_genes_dict)

#transposon insertion function
def transposon_insert(genomic_data,insert_positions,selected_insert_positions,insertions):
	rand_insert_positions = []
	insertion_sample=0
	pattern='TA'
	insert_sites=re.finditer(pattern,genomic_data)
	for insert in insert_sites:
		insert_position=insert.start()
		insert_positions.append(insert_position)
	while insertion_sample < insertions:
		insertion_sample+=1
		rand_insertion=random.choice(insert_positions)
		rand_insert_positions.append(rand_insertion)
	return(rand_insert_positions)	
	
#essential gene insert check function (3' insertion toleration)
def essential_gene_insert_check(essential_genes,rand_insert_positions):
	removed_positions=[]
	for gene in essential_genes:
		gene_location=essential_genes[gene]
		start=gene_location.start
		end=gene_location.end
		#three_prime_utr=int(end-100)
		for insert in rand_insert_positions:
			if insert>start and insert<end:
			#if insert>start and insert<three_prime_utr:
				rand_insert_positions.remove(insert)
	return(rand_insert_positions)	

#transposon downstream insert function
def transposon_seq_down(insert,genome,insert_lib,key):
	insert_start=int(insert+40)
	insert_end=int(insert+80)
	insert_positions = SeqFeature(FeatureLocation(insert_start,insert_end), type='domain')
	insert_seq=insert_positions.extract(genome)
	insert_lib[key]=insert_seq
	return(insert_lib)
	
#transposon upstream insert function
def transposon_seq_up(insert,genome,insert_lib,key):
	insert_start=int(insert+2)
	insert_end=int(insert+42)
	insert_positions = SeqFeature(FeatureLocation(insert_start,insert_end), type='domain')
	insert_seq=insert_positions.extract(genome)
	forward_seq=Seq(insert_seq)
	rev_complement=forward_seq.reverse_complement()
	insert_lib[key]=rev_complement
	return(insert_lib)
	
#fasta write function
def fasta_format(fasta_lib,seq,key,Mariner,seq_id):
	fasta_id='>Insert_Seq:'
	fasta_separator='|'
	seq_key=(fasta_id.strip() +str(key).strip()+fasta_separator.strip()+seq_id.strip())
	output_seq=(Mariner.strip()+seq.strip())
	fasta_lib[seq_key]=output_seq
	return(fasta_lib)
	
#write out fasta file
def fasta_file_write(fasta_lib,seq_id,insertions):
	outfile_1_name=('TraDIS_experiment_'.strip()+seq_id.strip()+'_'.strip()+str(insertions).strip()+'_transposons.fa'.strip())
	out_file_1=open(outfile_1_name,'w')
	for Id,seq in fasta_lib.items():
		out_file_1.write(Id)
		out_file_1.write('\n')
		out_file_1.write(str(seq))
		out_file_1.write('\n')
	out_file_1.close()	

	
#write out essential genes
def essential_gene_out(essential_genes,seq_id,insertions):
	outfile_2_name=('simulation_gene_locus_tags_'.strip()+seq_id.strip()+'_'.strip()+str(insertions).strip()+'_'.strip()+'.txt'.strip())
	out_file_2=open(outfile_2_name,'w')
	essential_genes=sorted(essential_genes)
	for locus_tag in essential_genes:
		out_file_2.write(locus_tag)
		out_file_2.write('\n')
	out_file_2.close()
