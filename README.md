# SV-genotyping
Evaluation of long-read SV genotyping methods

**1. Generation of simulated data**

```
#Download human reference genome
for i in {1..22} X Y;do echo $i;wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr${i}.fa.gz;done
zcat *fa.gz >> hg38_1_22XY.fa

#Simulated SV set
reference=hg38_1_22XY.fa
SURVIVOR simSV ${reference} parameter_file 0.01 0 sim_output

#Get variant sequences of INS and DEL
sed 's/>//g' simSV_output.insertions.fa|xargs -n 2|sed 's/\s/\t/g' > sim_output.insertions.txt
i=sim_output.insertions.txt
vcf=sim_output.vcf
awk 'ARGIND==1{a[$1]=$0}ARGIND>1&&($1 in a){print $0"\t"a[$1]}' $i \
<(perl -lane 'if(/#/){next}else{$CP=$F[0]."_$F[1]";print join("\t",$CP,@F[2..@F])}' $vcf)|perl -lane '{@CP=split /_/,$F[0];print join("\t",$CP[0],$CP[1],@F[1..2],$F[10],@F[4..8])}' - > sim_output_INS_SEQ.vcf
grep -E "#|SVTYPE=DEL" sim_output.vcf > sim_output_DEL.vcf
sniffles -n -1 -s 1 -m Sim_PB_CLR_minimap2.bam -v sim_output_DEL_SEQ.vcf --Ivcf sim_output_DEL.vcf
cat sim_output_INS_SEQ.vcf <(grep "SVTYPE=DEL" sim_output_DEL_SEQ.vcf) <(grep -E "SVTYPE=DUP|SVTYPE=INV|SVTYPE=TRA" sim_output.vcf) > sim_output_5SV.vcf

#Simulated PB CLR data
PB_error_profile=/SURVIVOR_v1.0.7/SURVIVOR/HG002_Pac_error_profile_bwa.txt
SURVIVOR simreads sim_output.fasta ${PB_error_profile} 30 PB_CLR_sim_read_30x.fasta
```

