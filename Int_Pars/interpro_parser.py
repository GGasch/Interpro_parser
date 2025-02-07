#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import argparse
import pandas as pd
import statistics


def get_args():

	parser=argparse.ArgumentParser(description ='Based on fasta and tsv file downloaded from Intepro, extract the portion of the sequence that matches the query domain')
	parser.add_argument('-f','--fasta', type=str,
				help='the fasta file with the full sequences as downloaded from Interpro')
	parser.add_argument('-t','--tsv',type=str,
				help='the tsv file with the boundaries of the hits downloaded from Interpro')
	parser.add_argument('-q','--query',type=str,
				help='the name of the InterPro family we are looking to extract')
	parser.add_argument('--trim', action="store_true",
				help='Remove entries that are likely to be truncated, by removing the sequence who match the IPR on less than AVG - 3SD of the IPR (Computed on the tsv)')
	args=parser.parse_args()
	
	return args



def stat (df) :
	#list_len = sum(df["Matches"].str.split(",").apply(lambda x: [int(i.split('..')[1])-int(i.split('..')[0]) for i in x]).tolist(),[])
	#Previous line if perfectly correct, but awfully nested. For the sake of futur reader, if any, I will denest it

	df_bound = df["Matches"].str.split(",")
	df_len=df_bound.apply(lambda x: [int(i.split('..')[1])-int(i.split('..')[0]) for i in x])
	list_len=sum(df_len.tolist(),[]) #Transfer the dataframe in a list. The sum is here to denest the list
	avg_len=statistics.mean(list_len)
	sd_len=statistics.stdev(list_len)
	return avg_len,sd_len



def seq_extract_no_trim (fasta_file,df,query) :

	list_record_crpd=[]

	with open (fasta_file) as fasta :
		for record in tqdm(SeqIO.parse(fasta,"fasta")):
			name=record.id[:record.id.find('|')]
			if name in df['Accession'].values : #Here to correct a small bug where Interpro gives fasta sequence that are absent from the tsv
				bounds_0= df.loc[df['Accession']==name, 'Matches'].iloc[0].split(',')#Extract the bounds (upper and lower) of all the domains present in the sequence
				bounds_1=[j.split('..') for j in bounds_0]#Turn the couple into a list of number that are easier to handle
				for bounds in bounds_1:
					seq_crpd=record.seq[int(bounds[0])-1:int(bounds[1])]#Crop the sequence base on the bounds
					record_crpd=SeqRecord(seq_crpd,record.id,"","")#Create the element of the future fasta file
					list_record_crpd.append(record_crpd)

	SeqIO.write(list_record_crpd,query+".fasta","fasta")




def seq_extract_trim(fasta_file, df, query,avg,sd):
	list_record_crpd = []

	with open(fasta_file) as fasta:
		for record in tqdm(SeqIO.parse(fasta, "fasta")):
			name = record.id[:record.id.find('|')]
			if name in df[
				'Accession'].values:  # Here to correct a small bug where Interpro gives fasta sequence that are absent from the tsv
				bounds_0 = df.loc[df['Accession'] == name, 'Matches'].iloc[0].split(',')  # Extract the bounds (upper and lower) of all the domains present in the sequence
				bounds_1 = [j.split('..') for j in bounds_0]  # Turn the couple into a list of number that are easier to handle
				for bounds in bounds_1 :
					if int(bounds[1])-int(bounds[0])>= avg-2*sd :
						seq_crpd = record.seq[int(bounds[0]) - 1:int(bounds[1])]  # Crop the sequence base on the bounds
						record_crpd = SeqRecord(seq_crpd, record.id, "", "")  # Create the element of the future fasta file
						list_record_crpd.append(record_crpd)

	SeqIO.write(list_record_crpd, query + ".fasta", "fasta")



def main():

	args=get_args()

	df_in = pd.read_csv(args.tsv, sep='\t')#Transfer the tsv file into a pandas file
	if args.trim :
		avg, sd = stat(df_in)
		print(avg, sd)
		seq_extract_trim(args.fasta, df_in, args.query, avg, sd)
	else:
		seq_extract_no_trim(args.fasta,df_in,args.query)








if __name__=='__main__':
	main()
