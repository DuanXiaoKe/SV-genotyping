#!/usr/bin/python3

import sys
import random as rd
#import numpy as np
from pyfaidx import Fasta
import fileinput
import re

# Load reference genome
fa = Fasta(sys.argv[1])
#ATCG = np.array(["A", "T", "C", "G"])
#error_seq = "".join(ATCG[[rd.randint(0,3) for i in range(size)]])
#error_num = rd.randint(1,3)
#error_seq = "".join(ATCG[[rd.randint(0,3) for i in range(error_num)]])

size=int(sys.argv[3])

out_file = open(sys.argv[4], 'w')

for line in fileinput.input(sys.argv[2]):
    # if header, print and got to next line
    if line[0] == '#':
        line = line.rstrip()
        print(line,sep="\t",file=out_file)
        continue
    # else parse variant record
    line = line.rstrip().split('\t')
    chrom = line[0]
    pos = int(line[1])
    ref = line[3]
    alt = line[4]
    svtype=re.search( r'SVTYPE=([A-Z]+)',line[7],re.I).group(1)
    shift_lr = rd.choice(["left","right"])
    if shift_lr == "left":
        if svtype =="DUP" or svtype =="INV":
            ref_seq_num = pos - size
            nuc = fa[str(chrom)][ref_seq_num -1]
            line[3] = nuc.seq
            #line[4] = nuc.seq + str(alt)[1:]
            line[1] = pos - size
            LEN=int(re.search( r'SVLEN=(\d+)',line[7],re.I).group(1))
            SVLEN = "SVLEN=" + str(LEN)
            SVTYPE = "SVTYPE=" + str(svtype)
            END = "END=" + str(pos + LEN - size)
            line[7] = ";".join([SVTYPE,SVLEN,END])
        else:
            ref_seq_num = pos - size
            nuc = fa[str(chrom)][ref_seq_num -1]
            line[3] = nuc.seq
            #line[4] = nuc.seq + str(alt)[1:]
            line[1] = pos - size
            chr2=re.search( r'CHR2=(\d+)',line[7],re.I).group(1)
            #SVLEN = "SVLEN=" + str(LEN)
            #SVTYPE = "SVTYPE=" + str(svtype)
            CHR2 = "CHR2=" + str(chr2)
            END = "END=" + str(pos - size)
            line[7] = ";".join(["SVTYPE=BND",END,CHR2])
    else:
        if svtype =="DUP" or svtype =="INV":
            ref_seq_num = pos + size
            nuc = fa[str(chrom)][ref_seq_num -1]
            line[3] = nuc.seq
            #line[4] = nuc.seq + str(alt)[1:]
            line[1] = pos + size
            LEN=int(re.search( r'SVLEN=(\d+)',line[7],re.I).group(1))
            SVLEN = "SVLEN=" + str(LEN)
            SVTYPE = "SVTYPE=" + str(svtype)
            END = "END=" + str(pos + LEN + size)
            line[7] = ";".join([SVTYPE,SVLEN,END])
        else:
            ref_seq_num = pos + size
            nuc = fa[str(chrom)][ref_seq_num -1]
            line[3] = nuc.seq
            #line[4] = nuc.seq + str(alt)[1:]
            line[1] = pos + size
            chr2=re.search( r'CHR2=(\d+)',line[7],re.I).group(1)
            #SVLEN = "SVLEN=" + str(LEN)
            #SVTYPE = "SVTYPE=" + str(svtype)
            CHR2 = "CHR2=" + str(chr2)
            END = "END=" + str(pos + size)
            line[7] = ";".join(["SVTYPE=BND",END,CHR2])
    print(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],sep="\t",file=out_file)
out_file.close()
print("Done......")
