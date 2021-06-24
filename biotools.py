#!/usr/bin/env python3

import sys
import gzip
import random

def read_fasta(filename):
	name = None
	seqs = []
	
	fp = None
	if filename == '-':
		fp = sys.stdin
	elif filename.endswith('.gz'):
		fp = gzip.open(filename, 'rt')
	else:
		fp = open(filename)

	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()

def anti(seq):                                                      # Function returns the complementary strand of DNA
    comp = ''
    for nt in seq:
        if nt == 'A': comp = 'T' + comp
        elif nt == 'C': comp = 'G' + comp
        elif nt == 'G': comp = 'C' + comp
        elif nt == 'T': comp = 'A' + comp
        else: comp = 'X' + comp
    return comp

def gc(seq):
	count = 0
	for nt in seq:
		if nt == 'G' or nt == 'C':
			count += 1
	return count / len(seq)
	
def randseq(l, gc):
	dna = []
	for i in range(l):
		r = random.random()
		if r < gc:
			r = random.random()
			if r < 0.5: dna.append('G')
			else:		dna.append('C')
		else:
			r = random.random()
			if r < 0.5: dna.append('A')
			else:		dna.append('T')
	return ''.join(dna)
	
def hydro(seq):
	kd = 0
	for aa in seq:
		if aa == 'I': kd += 4.5
		elif aa == 'V': kd += 4.2
		elif aa == 'L': kd += 3.8
		elif aa == 'F': kd += 2.8
		elif aa == 'C': kd += 2.5
		elif aa == 'M': kd += 1.9
		elif aa == 'A': kd += 1.8
		elif aa == 'G': kd += -0.4
		elif aa == 'T': kd += -0.7
		elif aa == 'S': kd += -0.8
		elif aa == 'W': kd += -0.9
		elif aa == 'Y': kd += -1.3
		elif aa == 'P': kd += -1.6
		elif aa == 'H': kd += -3.2
		elif aa == 'E': kd += -3.5
		elif aa == 'Q': kd += -3.5
		elif aa == 'D': kd += -3.5
		elif aa == 'N': kd += -3.5
		elif aa == 'K': kd += -3.9
		elif aa == 'R': kd += -4.5
	return kd / len(seq)
	
def skew(seq):                                                      #(g-c)/(g+c)
	g_count = 0
	c_count = 0
	for nt in seq:
		if nt == 'G': g_count += 1
		elif nt == 'C': c_count += 1
	return (g_count - c_count) / (g_count + c_count)

def find_orfs(seq, minprot):                                        # Identifying and storing all possible open reading frames
    atgs = []
    orfs = []
    stops = {}
    for i in range(len(seq) -2):
        if seq[i:i+3] == 'ATG':                                     # Genes begin with ATG
            atgs.append(i)
    for atg in atgs:
        stop = None
        for i in range(atg, len(seq) -2, 3):
            codon = seq[i:i+3]
            if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':  # Genes end with stop codon
                stop = i 
                if stop not in stops:
                    stops[stop] = True                              
                    orf = seq[atg:stop+3]
                    if len(orf)/3 > minprot:
                        orfs.append(orf)
                break
    return orfs
    
    
def longest_orf(seq):
	atgs = []
	for i in range(len(seq) -2):
		if seq[i:i+3] == 'ATG': atgs.append(i)
	max_len = 0
	max_seq = None
	for atg in atgs:
		stop = None
		for i in range(atg, len(seq) -2, 3):
			codon = seq[i:i+3]
			if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
				stop = i 
				break
		if stop != None:
			cds_len = stop - atg +3
			if cds_len > max_len:
				max_len = cds_len
				max_seq = seq[atg:atg+cds_len]
	return max_seq
	if max_seq == None: return None
	
def translate(seq):
	assert(len(seq) % 3 == 0)
	pro = []
	for i in range(0, len(seq), 3):
		codon = seq[i:i+3]
		if codon == 'AAA': pro.append('K')
		elif codon == 'AAC': pro.append('N')
		elif codon == 'AAT': pro.append('N')
		elif codon == 'AAG': pro.append('K')
		elif codon == 'ACA': pro.append('T')
		elif codon == 'ACC': pro.append('T')
		elif codon == 'ACT': pro.append('T')
		elif codon == 'ACG': pro.append('T')
		elif codon == 'AGA': pro.append('R')
		elif codon == 'AGC': pro.append('S')
		elif codon == 'AGT': pro.append('S')
		elif codon == 'AGG': pro.append('R')
		elif codon == 'ATA': pro.append('I')
		elif codon == 'ATC': pro.append('I')
		elif codon == 'ATT': pro.append('I')
		elif codon == 'ATG': pro.append('M')

		elif codon == 'CAC': pro.append('H')
		elif codon == 'CAT': pro.append('H')
		elif codon == 'CAG': pro.append('Q')
		elif codon == 'CAA': pro.append('Q')
		elif codon == 'CCA': pro.append('P')
		elif codon == 'CCC': pro.append('P')
		elif codon == 'CCT': pro.append('P')
		elif codon == 'CCG': pro.append('P')
		elif codon == 'CGA': pro.append('R')
		elif codon == 'CGC': pro.append('R')
		elif codon == 'CGT': pro.append('R')
		elif codon == 'CGG': pro.append('R')
		elif codon == 'CTA': pro.append('L')
		elif codon == 'CTC': pro.append('L')
		elif codon == 'CTT': pro.append('L')
		elif codon == 'CTG': pro.append('L')

		elif codon == 'GAC': pro.append('D')
		elif codon == 'GAA': pro.append('E')
		elif codon == 'GAT': pro.append('D')
		elif codon == 'GAG': pro.append('E')
		elif codon == 'GCA': pro.append('A')
		elif codon == 'GCC': pro.append('A')
		elif codon == 'GCT': pro.append('A')
		elif codon == 'GCG': pro.append('A')
		elif codon == 'GGA': pro.append('G')
		elif codon == 'GGC': pro.append('G')
		elif codon == 'GGT': pro.append('G')
		elif codon == 'GGG': pro.append('G')
		elif codon == 'GTA': pro.append('V')
		elif codon == 'GTC': pro.append('V')
		elif codon == 'GTT': pro.append('V')
		elif codon == 'GTG': pro.append('V')
		
		elif codon == 'TAC': pro.append('Y')
		elif codon == 'TAA': pro.append('X')
		elif codon == 'TAT': pro.append('Y')
		elif codon == 'TAG': pro.append('X')
		elif codon == 'TCA': pro.append('S')
		elif codon == 'TCC': pro.append('S')
		elif codon == 'TCT': pro.append('S')
		elif codon == 'TCG': pro.append('S')
		elif codon == 'TGA': pro.append('X')
		elif codon == 'TGC': pro.append('C')
		elif codon == 'TGT': pro.append('C')
		elif codon == 'TGG': pro.append('W')
		elif codon == 'TTA': pro.append('L')
		elif codon == 'TTC': pro.append('F')
		elif codon == 'TTT': pro.append('F')
		elif codon == 'TTG': pro.append('L')
		
		else: pro.append(X)
	return ''.join(pro)

def entropy(seq, w):
	for i in range(len(seq) -w +1):
		win = seq[i: i+ w]
		a_ct = 0.00001
		c_ct = 0.00001
		g_ct = 0.00001
		t_ct = 0.00001
		for nt in win:
			p = []
			if nt == 'A': a_ct += 1
			elif nt == 'C': c_ct += 1
			elif nt == 'G': g_ct += 1
			elif nt == 'T': t_ct += 1
			aprob = float(a_ct/w)
			p.append(aprob)
			cprob = float(c_ct/w)
			p.append(cprob)
			gprob = float(g_ct/w)
			p.append(gprob)
			tprob = float(t_ct/w)
			p.append(tprob)
			h = 0
		for i in range(len(p)):
			h -= p[i] * math.log2(p[i])
			if h < 0.0001:
				h = 0
	return h
	
aa = {
	'AAA' : 'K',	'AAC' : 'N',	'AAG' : 'K',	'AAT' : 'N',
	'ACA' : 'T',	'ACC' : 'T',	'ACG' : 'T',	'ACT' : 'T',
	'AGA' : 'R',	'AGC' : 'S',	'AGG' : 'R',	'AGT' : 'S',
	'ATA' : 'I',	'ATC' : 'I',	'ATG' : 'M',	'ATT' : 'I',
	'CAA' : 'Q',	'CAC' : 'H',	'CAG' : 'Q',	'CAT' : 'H',
	'CCA' : 'P',	'CCC' : 'P',	'CCG' : 'P',	'CCT' : 'P',
	'CGA' : 'R',	'CGC' : 'R',	'CGG' : 'R',	'CGT' : 'R',
	'CTA' : 'L',	'CTC' : 'L',	'CTG' : 'L',	'CTT' : 'L',
	'GAA' : 'E',	'GAC' : 'D',	'GAG' : 'E',	'GAT' : 'D',
	'GCA' : 'A',	'GCC' : 'A',	'GCG' : 'A',	'GCT' : 'A',
	'GGA' : 'G',	'GGC' : 'G',	'GGG' : 'G',	'GGT' : 'G',
	'GTA' : 'V',	'GTC' : 'V',	'GTG' : 'V',	'GTT' : 'V',
	'TAA' : '*',	'TAC' : 'Y',	'TAG' : '*',	'TAT' : 'Y',
	'TCA' : 'S',	'TCC' : 'S',	'TCG' : 'S',	'TCT' : 'S',
	'TGA' : '*',	'TGC' : 'C',	'TGG' : 'W',	'TGT' : 'C',
	'TTA' : 'L',	'TTC' : 'F',	'TTG' : 'L',	'TTT' : 'F',
}

kdscale = {
	'I' : 4.5,	'T' : -0.7,	'D' : -3.5,
	'V' : 4.2,	'S' : -0.8,	'N' : -3.5,
	'L' : 3.8,	'W' : -0.9,	'K' : -3.9,
	'F' : 2.8,	'Y' : -1.3,	'R' : -4.5,
	'C' : 2.5,	'P' : -1.6,'M' : 1.9,
	'H' : -3.2,	'A' : 1.8,	'E' : -3.5,
	'G' : -0.4,	'Q' : -3.5,
	}
ISscale = {
	'I' : -0.31,	'T' : 0.14,	'D' : 1.23,
	'V' : 0.07,	'S' : 0.13,	'N' : 0.42,
	'L' : -0.56,	'W' : -1.85,	'K' : 0.99,
	'F' : -1.13,	'Y' : -0.94,	'R' : 0.81,
	'C' : -0.24,	'P' : 0.45,	'M' : -0.23,
	'H' : 0.96,	'A' : 0.17,	'E' : 2.02,
	'G' : 0.01,	'Q' : 0.58,
	}
OSscale = { 
	'I' : -1.12,	'T' : 0.25,	'D' : 3.64,
	'V' : -0.46,	'S' : 0.46,	'N' : 0.85,
	'L' : -1.25,	'W' : -2.09,	'K' : 2.80,
	'F' : -1.71,	'Y' : -0.71,	'R' : 1.81,
	'C' : -0.02,	'P' : 0.14,	'M' : -0.67,
	'H' : 2.33,	'A' : 0.50,	'E' : 3.63,
	'G' : 1.15,	'Q' : 0.77,
	}
IS_OSscale = { 
	'I' : -0.81,	'T' : 0.11,	'D' : 2.41,
	'V' : -0.53,	'S' : 0.33,	'N' : 0.43,
	'L' : -0.69,	'W' : -0.24,	'K' : 1.81,
	'F' : -0.58,	'Y' : 0.23,	'R' : 1.00,
	'C' : 0.22,	'P' : -0.31,	'M' : -0.44,
	'H' : 1.37,	'A' : 0.33,	'E' : 1.61,
	'G' : 1.14,	'Q' : 0.19,
	}

def cal_hydrophobicity(seq, method, w):
	if method == 'KD': scale = kdscale
	elif method == 'IS': scale = ISscale
	elif method == 'OS': scale = OSscale
	elif method == 'IS+OS': scale = IS_OSscale
	else: 
		print("Method ", method, " not found.")
	score = 0
	for aa in seq:
		if aa in scale: score += scale[aa]
	return score