#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import argparse
import pandas as pd

parser=argparse.ArgumentParser(description ='Based on an interpro protein file, and the associated tsv file, extract the sequence in the protein matching the interpro profile')
parser.add_argument('-i','--input_file', type=str,
			help='the fasta file with the full sequences as downloaded from Interpro')
parser.add_argument('-t','--tsv_file',type=str,
			help='the tsv file with the boundaries of the hits downloaded from Interpro')
parser.add_argument('-n','--name',type=str,
			help='the name of the domain')
args=parser.parse_args()

df_in = pd.read_csv(args.tsv_file, sep='\t')
df_out=pd.DataFrame()
list_record_crpd=[]

with open (args.input_file) as fasta :
	for record in tqdm(SeqIO.parse(fasta,"fasta")):
		name=record.id[:record.id.find('|')]
		bounds_0= df_in.loc[df_in['Accession']==name, 'Matches'].iloc[0].split(',')
		bounds_1=[j.split('..') for j in bounds_0]
		for bounds in bounds_1:
			seq_crpd=record.seq[int(bounds[0])-1:int(bounds[1])]
			record_crpd=SeqRecord(seq_crpd,record.id,"","")
			list_record_crpd.append(record_crpd)
SeqIO.write(list_record_crpd,args.name+".fasta","fasta")
