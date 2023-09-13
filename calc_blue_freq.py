#!/usr/bin/env python3 
import numpy
import gzip
import copy
import argparse

parser = argparse.ArgumentParser(prog="calc_blue_freq.py",description='Calculate BLUE of allele frequency. Can handle missing data. Outputs as a GARLIC .freq file.')

parser.add_argument('-v', '--vcf', type=str, required=True, metavar='VCF_FILE', help='Gzipped VCF file')
parser.add_argument('-k', '--king', type=str, required=True, metavar='KIN_MAT', help='Square kinship matrix. e.g. from plink --make-king square')
parser.add_argument('-i', '--ids', type=str, required=True, metavar='KING_IDS', help='List of individual IDs in order of kinship matrix. Should also be same order of samples in VCF. Assumes first row is a header.')

args = parser.parse_args()

#reads in kinship matrix and id list direct output of plink '--make-king square'
def kin_matrix(idfile, kingfile):
	ids = []
	with open(idfile) as f:
		for line in f:
			ids.append(line.split('\n')[0])
	ids.pop(0) #get rid if first row, it is a label
	kin = numpy.empty((len(ids), len(ids)))
	with open(kingfile) as f:
		counter=0
		for line in f:
			#set negative entries to 0, not valid kin coefs, related to population stucture in sample
			temp = [float(x.replace('\n', '')) if float(x.replace('\n', '')) > 0 else 0.0 for x in line.split('\t')]
			kin[counter]=temp
			counter=counter+1
	return ids, kin

#kin: kinship matrix (symmetric non-negative invertable)
#ids: list of individual ids (len(ids) == dim of kin)
#to_delete: list of indices to delete, non-repeating, assumed to be in descending order (e.g., [4,3,0] not [0,3,4])
#all_weights: dict storing weights for a given submatrix (keyed by concatenation of ids, sep = "."")
#all_weights is assumed to be seeded with the weights from the full matrix, this could be changed
def get_weights(kin,ids,to_delete,all_weights):
	#nothing to delete, return all_weights of everyone
	if len(to_delete) == 0:
		return all_weights[".".join(ids)], all_weights
	#list to hold new ids
	newids = copy.deepcopy(ids)
	#this is why to_delete is assumed to be in descending order
	#if not, could go OOB or drop the wrong individual
	for i in to_delete:
		newids.pop(i)
	#create key for the subset
	label = ".".join(newids)
	#did we already compute it?
	if label in all_weights:
		return all_weights[label], all_weights
	#guess we did not
	else:
		#take kin and delete rows and columns
		kin0 = numpy.delete(kin,to_delete,axis=0)
		kin1 = numpy.delete(kin0,to_delete,axis=1)
		#calculate weights McPeek et al 2004 Biometrics
		kin_inv = numpy.linalg.inv(kin1)
		ones_vector = numpy.ones((kin1.shape[0], 1))
		num_sum = kin_inv.sum(axis=0)
		denominator = ones_vector.T @ kin_inv @ ones_vector
		w = num_sum/denominator
		all_weights[label] = w
		return w, all_weights


idfile=args.ids #"dingo_NGSD_20samples_gatk_merged_chraut_variants_biallelic_snps_filtered_nomono_het.king.id"
kingfile=args.king #"dingo_NGSD_20samples_gatk_merged_chraut_variants_biallelic_snps_filtered_nomono_het.king"
vcffile=args.vcf #"dingo_NGSD_20samples_gatk_merged_chraut_variants_biallelic_snps_filtered_nomono_het.vcf.gz"

ids, kin = kin_matrix(idfile, kingfile)

kin_inv = numpy.linalg.inv(kin)
ones_vector = numpy.ones((kin.shape[0], 1))
num_sum = kin_inv.sum(axis=0)
denominator = ones_vector.T @ kin_inv @ ones_vector
w = num_sum/denominator

all_weights = {}
label=".".join(ids)
all_weights[label] = w

print("CHR SNP POS ALLELE FREQ")
with gzip.open(vcffile,'rt') as fin:
	for line in fin:
		if line[0] == '#':
			continue
		fields = line.strip().split()
		info=fields[0:9]
		gts=fields[9:]
		counts = []
		to_delete = []
		for i, gt in enumerate(gts):
			gt,fmt = gt.split(':',1)
			if gt == '.' or gt == './.' or gt == '.|.':
				to_delete.append(i)
				continue
			else:
				gt = int(gt[0])+int(gt[2])
			counts.append(gt)
		if len(to_delete) > 0:
			to_delete.reverse()
		npcounts = numpy.asarray(counts)/2
		w, all_weights = get_weights(kin,ids,to_delete,all_weights)
		npcounts = numpy.sum(w*npcounts)

		print(info[0],end=" ")
		print(info[2],end=" ")
		print(info[1],end=" ")
		print(info[4],end=" ")
		print(npcounts)		

