import numpy as np
import pandas as pd

# Nussinov algorithm for RNA secondary structure prediction
# predict the secondary structure of an RNA sequence and then visualize the predicted structure in dot-bracket notation


# check if RNA nucleotides are Watson-Crick base pairs
def couple(pair):
    pairs = {"A": "U", "U": "A", "C": "G", "G": "C"}
    if pair in pairs.items():
        return True
    
    return False

# fill the matrix accroding to the Nussinov algorithm
# matrix used to find the optimal secondary structure

def fill(nm, rna):
    minLoopLength = 0;
    for k in range(1, len(rna)):
        for i in range(len(rna) - k):
            j = i + k  
            if j - i >= minLoopLength:
                down = nm[i + 1][j] # 1st rule
                left = nm[i][j - 1] # 2nd rule
                diag = nm[i + 1][j - 1] + couple((rna[i], rna[j])) # 3rd rule

                rc = max([nm[i][t] + nm[t + 1][j] for t in range(i, j)]) # 4th rule
                nm[i][j] = max(down, left, diag, rc)
            else:
                nm[i][j] = 0
    return nm

# traceback through complete Nussinov matrix to find optimial RNA secondary structure solution through maximum base-pairs
def traceback(nm, rna, fold, i, L):
	j = L
	if i < j:
		if nm[i][j] == nm[i + 1][j]: # 1st rule
			traceback(nm, rna, fold, i + 1, j)
		elif nm[i][j] == nm[i][j - 1]: # 2nd rule
			traceback(nm, rna, fold, i, j - 1)
		elif nm[i][j] == nm[i + 1][j - 1] + couple((rna[i], rna[j])): # 3rd rule
			fold.append((i, j))
			traceback(nm, rna, fold, i + 1, j - 1)
		else:
			for k in range(i + 1, j - 1):
				if nm[i][j] == nm[i, k] + nm[k + 1][j]: # 4th rule
					traceback(nm, rna, fold, i, k)
					traceback(nm, rna, fold, k + 1, j)
					break
				
	return fold

# generate the dot-bracket notation 

def dotWrite(rna, fold):
	dot = ["." for i in range(len(rna))]

	for s in fold:
		#print(min(s), max(s))
		dot[min(s)] = "("
		dot[max(s)] = ")"

	return "".join(dot)


def matrixInitialisation(rna):
	M = len(rna)

	# init matrix
	nm = np.empty([M, M]) # np.empty to save space instead of np.zeroes
	nm[:] = np.NAN

	# init diagonals to 0
	nm[range(M), range(M)] = 0
	nm[range(1, len(rna)), range(len(rna) - 1)] = 0

	return nm


rna = """GCGGCAACAGCGGGGCCGATGTGTAGTTGGTGACTGCCTCTCCAGATGCTGAGGTGCCTG
TATCATTGGCACAGGCCAGTGCTGAACCGTAGGTGGAGTAGGCTGTGCCTTCTGAAGCAG
TATCTATTCACAATGAAGTTGCAGTCTCCCGAATTCCAGTCACTTTTCACAGAAGGACTG
AAGAGTCTGACAGAATTATTTGTCAAAGAGAATCACGAATTAAGAATAGCAGGAGGAGCA
GTGAGGGATTTATTAAATGGAGTAAAGCCTCAGGATATAGATTTTGCCACCACTGCTACC
CCTACTCAAATGAAGGAGATGTTTCAGTCGGCTGGGATTCGGATGATAAACAACAGAGGA
GAAAAGCACGGAACAATTACTGCCAGGCTTCATGAAGAAAATTTTGAGATTACTACACTA
CGGATTGATGTCACCACTGATGGAAGACATGCTGAGGTAGAATTTACAACTGACTGGCAG
AAAGATGCGGAACGCAGAGATCTCACTATAAATTCTATGTTTTTAGGTTTTGATGGCACT
TTATTTGACTACTTTAATGGTTATGAAGATTTAAAAAATAAGAAAGTTAGATTTGTTGGA
CATGCTAAACAGAGAATACAAGAGGATTATCTTAGAATTTTAAGATACTTCAGGTTTTAT
GGGAGAATTGTAGACAAACCTGGTGACCATGATCCTGAGACTTTGGAAGCAATTGCAGAA
AATGCAAAAGGCTTGGCTGGAATATCAGGAGAAAGGATTTGGGTGGAACTGAAAAAAATT
CTTGTTGGTAACCATGTAAATCATTTGATTCACCTTATCTATGATCTTGATGTGGCTCCT
TATATAGGTTTACCTGCTAATGCAAGTTTAGAAGAATTTGACAAAGTCAGTAAAAATGTT
GATGGTTTTTCACCAAAGCCAGTGACTCTTTTGGCCTCATTATTCAAAGTACAAGATGAT
GTCACAAAATTGGATTTGAGGTTGAAGATCGCGAAAGAGGAGAAAAACCTTGGCTTATTT
ATAGTTAAAAATAGGAAAGATTTAATTAAAGCAACAGATAGTTCAGACCCATTGAAACCC
TATCAAGACTTCATTATAGATTCTAGGGAACCTGATGCAACTACTCGTGTATGTGAACTA
CTGAAGTACCAAGGAGAGCACTGTCTCCTAAAGGAAATGCAGCAGTGGTCCATTCCTCCA
TTTCCTGTAAGTGGCCATGACATCAGAAAAGTGGGCATTTCTTCAGGAAAAGAAATTGGG
GCTCTATTACAACAGTTGCGAGAACAGTGGAAAAAAAGTGGTTACCAAATGGAAAAAGAT
GAACTTCTGAGTTACATAAAGAAGACCTAAAACTGATGGCTACTAAAAAGCAGAGCATTT"""

rna = rna.replace("\n","")
nm = matrixInitialisation(rna)
nm = fill(nm, rna)
fold = []
sec = traceback(nm, rna, fold, 0, len(rna) - 1)
res = dotWrite(rna, fold)


# create a Pandas DataFrame df to display the Nussinov matrix,
# using the nucleotide names as row and column labels
names = [_ for _ in rna]
df = pd.DataFrame(nm, index=names, columns=names)

# print the Nussinov matrix, the original RNA sequence,
# and the predicted secondary structure in dot-bracket notation
print(df, "\n", rna, res)