#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import argparse
import pandas as pd


def get_args():

	parser=argparse.ArgumentParser(description ='Based on fasta and tsv file downloaded from Intepro, extract the portion of the sequence that matches the query domain')
	parser.add_argument('-f','--fasta', type=str,
				help='the fasta file with the full sequences as downloaded from Interpro')
	parser.add_argument('-t','--tsv',type=str,
				help='the tsv file with the boundaries of the hits downloaded from Interpro')
	parser.add_argument('-q','--query',type=str,
				help='the name of the InterPro family we are looking to extract')
	args=parser.parse_args()
	
	return args

def main():
	args=get_args()

	df_in = pd.read_csv(args.tsv, sep='\t')#Transfer the tsv file into a pandas file
	list_record_crpd=[]

	with open (args.fasta) as fasta :
		for record in tqdm(SeqIO.parse(fasta,"fasta")):
			name=record.id[:record.id.find('|')]
			if name in df_in['Accession'].values : #Here to correct a small bug where Interpro gives fasta sequence that are absent from the tsv
				bounds_0= df_in.loc[df_in['Accession']==name, 'Matches'].iloc[0].split(',')#Extract the bounds (upper and lower) of all the domains present in the sequence
				bounds_1=[j.split('..') for j in bounds_0]#Turn the couple into a list of number that are easier to handle
				for bounds in bounds_1:
					seq_crpd=record.seq[int(bounds[0])-1:int(bounds[1])]#Crop the sequence base on the bounds
					record_crpd=SeqRecord(seq_crpd,record.id,"","")#Create the element of the future fasta file
					list_record_crpd.append(record_crpd)
			else :
				with open("error.txt", "a") as file:
					file.write(f'{name}\n')
	SeqIO.write(list_record_crpd,args.query+".fasta","fasta")



if __name__=='__main__':
	main()
