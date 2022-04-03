#!/usr/bin/python3

import sys
import random as rd
from pyfaidx import Fasta
import fileinput
import re

fa = Fasta(sys.argv[1])
size=int(sys.argv[3])
out_file = open(sys.argv[4], 'w')

for line in fileinput.input(sys.argv[2]):
    if line[0] == '#':
        line = line.rstrip()
        print(line,sep="\t",file=out_file)
        continue
    line = line.rstrip().split('\t')
    chrom = line[0]
    pos = int(line[1])
    ref = line[3]
    alt = line[4]
    svtype=re.search( r'SVTYPE=([A-Z]+)',line[7],re.I).group(1)
    shift_lr = rd.choice(["left","right"])
    if shift_lr == "left":
        if svtype =="INS":
            ref_seq_num = pos - size
            nuc = fa[str(chrom)][ref_seq_num -1]
            line[3] = nuc.seq
            line[4] = nuc.seq + str(alt)[1:]
            line[1] = pos - size
            SVLEN = "SVLEN=" + str(len(line[4])-1)
            END = "END=" + str(pos - size)
            line[7] = ";".join(["SVTYPE=INS",SVLEN,END])
        else:
            ref_seq_num = pos - size
            ref_seq_num_1 = pos + len(line[3]) - 1 - size
            nuc = fa[str(chrom)][ref_seq_num:ref_seq_num_1]
            nuc_1 = fa[str(chrom)][ref_seq_num - 1]
            line[3] = nuc_1.seq + nuc.seq
            line[4] = nuc_1.seq
            line[1] = ref_seq_num
            SVLEN = "SVLEN=-" + str(len(line[3])-1)
            END = "END=" + str(pos - size + len(line[3])-1)
            line[7] = ";".join(["SVTYPE=DEL",SVLEN,END])
    else:
        if svtype =="INS":
            ref_seq_num = pos + size
            nuc = fa[str(chrom)][ref_seq_num -1]
            line[3] = nuc.seq
            line[4] = nuc.seq + str(alt)[1:]
            line[1] = pos + size
            SVLEN = "SVLEN=" + str(len(line[4])-1)
            END = "END=" + str(pos + size)
            line[7] = ";".join(["SVTYPE=INS",SVLEN,END])
        else:
            ref_seq_num = pos + size
            ref_seq_num_1 = pos + len(line[3]) - 1 + size
            nuc = fa[str(chrom)][ref_seq_num:ref_seq_num_1]
            nuc_1 = fa[str(chrom)][ref_seq_num - 1]
            line[3] = nuc_1.seq + nuc.seq
            line[4] = nuc_1.seq
            line[1] = ref_seq_num
            SVLEN = "SVLEN=-" + str(len(line[3])-1)
            END = "END=" + str(pos + size + len(line[3])-1)
            line[7] = ";".join(["SVTYPE=DEL",SVLEN,END])
    print(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],sep="\t",file=out_file)
out_file.close()
print("Done......")
