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
        print(line)
        continue
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
            line[1] = pos - size
            chr2=re.search( r'CHR2=([\w\d]+)',line[7],re.I).group(1)
            end=re.search( r'END=([\w\d]+)',line[7]).group(1)
            CHR2 = "CHR2=" + str(chr2)
            END = "END=" + str(int(end) - size)
            AEND = str(int(end) - size)
            line[7] = ";".join(["SVTYPE=BND;SVLEN=1",CHR2,END])
            ALT = re.sub(r':\d+', ":"+AEND, line[4])
            line[4] = ALT
    else:
        if svtype =="DUP" or svtype =="INV":
            ref_seq_num = pos + size
            nuc = fa[str(chrom)][ref_seq_num -1]
            line[3] = nuc.seq
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
            line[1] = pos + size
            chr2=re.search( r'CHR2=([\w\d]+)',line[7],re.I).group(1)
            end=re.search( r'END=([\w\d]+)',line[7]).group(1)
            CHR2 = "CHR2=" + str(chr2)
            END = "END=" + str(int(end) + size)
            AEND = str(int(end) + size)
            line[7] = ";".join(["SVTYPE=BND;SVLEN=1",CHR2,END])
            ALT = re.sub(r':\d+', ":"+AEND, line[4])
            line[4] = ALT
    print(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],sep="\t",file=out_file)
out_file.close()
print("Done......")
