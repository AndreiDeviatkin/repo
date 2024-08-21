import argparse
import os
import re
#from collections import defaultdict

from Bio import SeqIO

def yn2table(inputfile):
    outfile_name = os.path.splitext(inputfile)[0] + 'YN.txt'

    with open(inputfile) as file:
        t = 0 # flag that means reaching section with YN results
        lines_out = []
        ids_num = {}
        count = 0
        for line in file:
            if line == "\n":
                continue
            if line == "Nei & Gojobori 1986. dN/dS (dN, dS)\n":
                t = 'NG'
                continue
            if line == "(B) Yang & Nielsen (2000) method\n":
                t = ""
                continue
            if line == "seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE\n":
                t = "YN"
                line = line.strip(" ").replace("+-","")
                line = ' '.join(line.split()).replace(' ', '\t') + '\n'
                lines_out.append(line)
                continue

            if line == "(C) LWL85, LPB93 & LWLm methods\n":
                t="LM"
                continue
            if t == "NG":
                if  line != "(Note: This matrix is not used in later ML. analysis.\n" and line!="Use runmode = -2 for ML pairwise comparison.)\n":
                    count +=1
                    ids_num[str(count)] = line.split(' ')[0]
            if t == "YN":
                line = line.strip(" ").replace("+-","")
                line = ' '.join(line.split()).replace(' ', '\t')  + '\n'
                line_l = line.split('\t')
                num1 = line_l[0]
                num2 = line_l[1]
                line = '\t'.join([ids_num[num1], ids_num[num2]]+line_l[2:])
                lines_out.append(line)

    #print(ids_num)
    print("Number of ids {}".format(len(ids_num)))
    file.close()
    
    with open(outfile_name, 'w') as file:
        file.writelines(lines_out)
    file.close()


def change_yn_name(input_yn, input_al, type):

    outfile_name = os.path.splitext(input_yn)[0] + '_ed.txt'
    
    dict_names = {}
    with open(input_al) as file:
        seqs = list(SeqIO.parse(file, "fasta"))
        for seq in seqs:
            dict_names[seq.id.split('_')[-1]] = seq.id
    file.close()
    
    
    if type == 'yn00':
        sep= '\t'
    else:
        sep = ','
    
    lines_out = []
    with open(input_yn) as file:
        #print(file.readline())
        lines = file.readlines()
        if type == 'yn00':
            line0 = ','.join(['name1', 'name2'] + lines[0].split('\t')[2:])
        else:
            line0 = lines[0]
        lines_out.append(line0)
        for line in lines[1:]:
            line_l = line.split(",")
            line = ','.join([dict_names[line_l[0]], dict_names[line_l[1]]]+line_l[2:])
            lines_out.append(line)
    file.close()
    
    with open(outfile_name, 'w') as file:
        file.writelines(lines_out)
    file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file with yn00.exe", required=True)
    parser.add_argument("-al", "--al", type=str,
                        help="Alignment with full sequence names", required=False)
    args = parser.parse_args()

    yn2table(args.input_file)
    file_in = os.path.splitext(args.input_file)[0] + 'YN.txt'

    if args.al != None:
        change_yn_name(file_in, args.al,args.type)