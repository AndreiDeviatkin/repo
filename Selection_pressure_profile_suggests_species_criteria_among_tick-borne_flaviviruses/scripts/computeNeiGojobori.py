import argparse
from Bio import SeqIO
from Bio import AlignIO, Align
import os.path
'''
computeNeiGojobori.py takes alignment of nucleotide sequences in fasta format and outputs a table in csv format with
 *numbers of synonymous and non-synonymous nucleotide substitutions
 *numbers of  potentially synonymous and potentially nonsynonymous sites
 *p-distances (count of the number of synonymous differences (Sd) normalized using the possible number of synonymous sites (S)) 
 for each sequence pair using Nei-Gojobori method (Nei and Gojobori 1986).
 
 Two files with precomputed observed and potential synonymous and nonsynonymous differences between codons are required (codonPairs_SdNd.csv and codonPairs_SN.csv respectively).
 
'''
def calc_differences(rec1,rec2):
    '''
    Computes number of synonymous and non-synonymous nucleotide substitutions and the numbers of 
    potentially synonymous and potentially nonsynonymous sites in two sequences using Nei-Gojobori method 
    (Nei and Gojobori 1986).

    S_obs, N_obs - counts of the number of synonymous (Sd) and nonsynonymous (Nd) differences.

    S, N - possible number of synonymous sites and nonsynonymous sited (average number for pair)

    rec1,rec2 - SeqRecords to be compared
    '''
    c = 3 # codon length
    S_obs = 0
    N_obs = 0
    S = 0
    N = 0
    
    seq1 = str(rec1.seq).upper()
    seq2 = str(rec2.seq).upper()
    
    for i in range(0,len(rec1.seq),c):
        codon_s1 = seq1[i:i+c]
        codon_s2 = seq2[i:i+c]
        if len(codon_s1)%3 != 0:
            continue
        if '-' in codon_s1 or '-' in codon_s2:
            continue
        #print(codon_s1,codon_s2,SdNd[(codon_s1,codon_s2)])
        S_obs += SdNd[(codon_s1,codon_s2)]['S']
        N_obs += SdNd[(codon_s1,codon_s2)]['N']
        S +=SN[(codon_s1,codon_s2)]['S']
        N += SN[(codon_s1,codon_s2)]['N']
    return [rec1.id, rec2.id, S_obs, N_obs, S, N, S_obs/S, N_obs/N]

def precalculated_SN():
    ''' loads files with numbers of substitutions and potential synonymous and nonsynonymous substitutions for codon pairs
    '''
    SdNd = {}

    with open('codonPairs_SdNd.csv') as file:
        for line in file:
            line_l = line.strip('\n').split(',')
            SdNd[tuple((line_l[0],line_l[1]))] = {'S':float(line_l[2]),'N':float(line_l[3])}
            
    SN = {}

    with open('codonPairs_SN.csv') as file:
        for line in file:
            line_l = line.strip('\n').split(',')
            SN[tuple((line_l[0],line_l[1]))] = {'S':float(line_l[2]),'N':float(line_l[3])}
    return SdNd,SN


def complete_deletion(aln_object):
    '''
    Codon-wise deletion of columns with gaps in alignment object
    '''
    aln_to_analyse = Align.MultipleSeqAlignment([])
    length = len(aln_object[0])
    #length of codon
    c=3
    # iterates over codons in alignment
    for i in range(0, length, c):
        skip=0
        codon_slice = aln_object[:,i:i+c]
        # check presence of gap in each position of codon
        for pos in range(c):
            if '-' in codon_slice[:,pos]:
                skip=1
                break
        if skip == 0:
            if len(aln_to_analyse) != 0:
                aln_to_analyse = aln_to_analyse + codon_slice
            else:
                aln_to_analyse = codon_slice
        else:
            continue
    return aln_to_analyse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file with alignment in fasta format", required=True)
    parser.add_argument("-complete_del", "--complete_del", type=int,
                        help="Complete deletion or pairwise deletion. If complete deletion, 1.")
    parser.add_argument("-output", "--output_name", type=str,
                        help="Name of output table")
    args = parser.parse_args()

    # load precomputed numbers of synonymous and non-synonymous nucleotide substitutions and the numbers of 
    # potentially synonymous and potentially nonsynonymous sites for codon pairs
    SdNd,SN = precalculated_SN()
    
    # load alignment of sequences
    aln = AlignIO.read(args.input_file, 'fasta')
    if args.complete_del:
        print('Complete deletion')
        aln_to_analyse = complete_deletion(aln)
    else:
        aln_to_analyse = aln

    # strings of table with distances
    table_str = ['name1,name2,Sd,Nd,S,N,pS,pN\n']
    for i in range(len(aln_to_analyse)):
        for j in range(i+1,len(aln_to_analyse)):
            d = calc_differences(aln_to_analyse[i],aln_to_analyse[j])
            table_str.append(','.join(str(x) for x in d)+'\n')
    if not args.output_name:
        args.output_name = os.path.splitext(args.input_file)[0] + '_NeiGoj.csv'
    
    with open(args.output_name,'w') as file:
        file.writelines(table_str)
    file.close()
