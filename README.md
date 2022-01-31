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

**2. Acquisition of GIAB data**

(1) PB CLR data for HG002, HG003, and HG004

HG002

```
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/PacBio_fasta/MD5.txt
d=ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/PacBio_fasta
for i in `cat MD5.txt|awk '{print $2}'`;do wget -c $d/${i};done
zcat m*subreads.fasta.gz >> HG002_CLR.fasta
```

HG003

```
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/PacBio_MtSinai_NIST/PacBio_fasta/MD5.txt
d=ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/PacBio_MtSinai_NIST/PacBio_fasta
for i in `cat MD5.txt|awk '{print $2}'`;do wget -c $d/${i};done
zcat m*subreads.fasta.gz >> HG003_CLR.fasta
```

HG004

```
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/PacBio_MtSinai_NIST/PacBio_fasta/MD5.txt
d=ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/PacBio_MtSinai_NIST/PacBio_fasta
for i in `cat MD5.txt|awk '{print $2}'`;do wget -c $d/${i};done
zcat m*subreads.fasta.gz >> HG004_CLR.fasta
```

(2) PB CCS data for HG002

Copy the file names of HG002 PB CCS sequence data (available at ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/) to seq_fq_id. txt

```
data=seq_fq_id.txt
d=ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb
for i in `cut -f 1 $data`;do wget -c $d/$i;done
cat *.seq_data >> HG002_CCS.seq_data
```

(3) ONT data for HG002

```
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_1.seq_data.gz
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_2.seq_data.gz
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_3.seq_data.gz
zcat GM24385_*.seq_data.gz >> HG002_ONT.seq_data
```

(4) PB CLR data for HG005, HG006, and HG007

HG005

```
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/MtSinai_PacBio/PacBio_fasta/MD5.txt
d=ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/MtSinai_PacBio/PacBio_fasta
for i in `cat MD5.txt|awk '{print $2}'`;do wget -c $d/${i};done
zcat m*subreads.fasta.gz >> HG005_CLR.fasta
```

HG006

```
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG006_NA24694-huCA017E_father/PacBio_MtSinai/PacBio_fasta/MD5.txt
d=ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG006_NA24694-huCA017E_father/PacBio_MtSinai/PacBio_fasta
for i in `cat MD5.txt|awk '{print $2}'`;do wget -c $d/${i};done
zcat m*subreads.fasta.gz >> HG006_CLR.fasta*
```

HG007

```
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG007_NA24695-hu38168_mother/PacBio_MtSinai/PacBio_fasta/MD5.txt
d=ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG007_NA24695-hu38168_mother/PacBio_MtSinai/PacBio_fasta
for i in `cat MD5.txt|awk '{print $2}'`;do wget -c $d/${i};done
zcat m*subreads.fasta.gz >> HG007_CLR.fasta
```

(5) PB CCS data for HG005

Obtain SRA accession number (starting SRR) of HG005 PB CCS sequence data (available at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA540706) and save them in the HG005_CCS_SRR.txt

```
prefetch --output-directory . --option-file HG005_CCS_SRR.txt
for i in `ls *sra`;do fasterq-dump $i;done
cat *seq_data > HG005_CCS.seq_data
```

(6) HG002 SV set

```
wget -c https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
#HG002 Tier 1
awk -v OFS="\t" '{if(/#/){print}else if($7=="PASS"){print}}' <(zcat HG002_SVs_Tier1_v0.6.vcf.gz) > HG002_Tier1.vcf
#HG002 Tier 2
awk -v OFS="\t" '{if(/#/){print}else if($7=="ClusteredCalls"){print}}' <(zcat HG002_SVs_Tier1_v0.6.vcf.gz)|grep -E -v "\./" > HG002_Tier2.vcf
```

(7) HG005 SV set

```
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
wget -c https://github.com/PacificBiosciences/pbsv/blob/master/annotations/human_hs37d5.trf.bed
pbmm2 align -j 4 -J 4 --sort --preset CCS --sample HG005 hs37d5.fa.gz HG005_CCS.seq_data HG005_pbmm2.bam
pbsv discover --tandem-repeats human_hs37d5.trf.bed HG005_pbmm2.bam HG005_pbmm2.svsig.gz
pbsv call --ccs -A 3 -O 3 -P 20 --gt-min-reads 3 -j 4 -t DEL,INS,INV,DUP,BND hs37d5.fa.gz HG005_pbmm2.svsig.gz HG005_pbmm2_pbsv.vcf

#Employ the SVs that are > 50 bp, with the FILTER “PASS” tag, and determined genotype for the further analyses
cat <(grep "#" HG005_pbmm2_pbsv.vcf) <(perl -lane 'if(/SVTYPE=INV/){($END)=$F[7]=~m/;END=([^;]+)/;$LEN=$END-$F[1];$INFO=$F[7].";SVLEN=$LEN";print join("\t",@F[0..6],$INFO,@F[8..@F])}else{print}' $vcf|bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\tSVTYPE=%SVTYPE;SVLEN=%SVLEN;END=%END\tGT\t[%GT]\n'|perl -lane 'if(/SVTYPE=BND/){$F[7]=~s/SVLEN=.*?;//g;($chr2)=$F[4]=~m/([0-9A-Z]{1,2}):/;$CHR2="CHR2=".$chr2;print join("\t",@F[0..6],join(";",$F[7],$CHR2),@F[8..@F])}else{if(/SVLEN=(.*?);/){$LEN=abs($1)}if($LEN>=50){print}}') |grep -E -v "\.\/\."|grep -E -v "NC|GL|X:|Y:"|perl -lane 'if(/^[^#|\d]/){next}else{print}'|grep -E "#|PASS" > HG005_pbmm2_pbsv_5SV.vcf
```

**3.  Read alignment**

(1) Alignment

```
reference=hs37d5.fa.gz
seq_data=fasta or fastq file for simulated, HG002, or HG005
thread=4
sample=simulated, HG002 or HG005
map=minimap2 or ngmlr
```

```
#minimap2 alignment
#CLR:
minimap2 -ax map-pb ${reference} ${seq_data} --MD -Y -t ${thread} -R '@RG\tID:${sample}' -o ${sample}_${map}.sam
#CCS:
minimap2 -ax asm20 ${reference} ${seq_data} --MD -Y -t ${thread} -R '@RG\tID:${sample}' -o ${sample}_${map}.sam
#ONT:
minimap2 -ax map-ont ${reference} ${seq_data} -z 600,200 --MD -Y -t ${thread} -R '@RG\tID:${sample}' -o ${sample}_${map}.sam
```

```
#ngmlr alignment
#CLR/CCS:
ngmlr -t ${thread} -r ${reference} -q ${seq_data} -x pacbio -o ${sample}_${map}.sam
#ONT:
ngmlr -t ${thread} -r ${reference} -q ${seq_data} -x ont -o ${sample}_${map}.sam
```

(2) Alignment sorting and down-sampling

```
samtools view -buS -@ ${thread} ${sample}_${map}.sam | \
samtools sort -@ ${thread} -O BAM - > ${sample}_${map}.bam && samtools index ${sample}_${map}.bam

samtools view -bS -s ${ratio} ${sample}_${map}.bam > ${sample}_${map}_ratio.bam && samtools index ${sample}_${map}_${cov}.bam
HG002 PB CLR data: $ratio are 0.87(60x), 0.73(50x), 0.58(40x), 0.44(30x), 0.30(20x), 0.15(10x), and 0.075(5x), respectively.
HG005 PB CLR data: $ratio are 0.88(50x), 0.71(40x), 0.53(30x), 0.36(20x), 0.18(10x), and 0.09(5x), respectively.
```

**4. SV genotyping**

```
reference=hs37d5.fa.gz
sample=simulated, HG002, or HG005
platform=CLR, CCS, or ONT
map=minimap2 or ngmlr
cov=5x, 10x, 20x, 30x, 40x, 50x, or 60x
bam=${sample}_${platform}_${map}_${cov}.bam
input_vcf=sim_output_5SV.vcf, HG002_Tier1.vcf, HG002_Tier2.vcf, or HG005_pbmm2_pbsv_5SV.vcf
thread=4
```

(1) cuteSV

```
#CLR
cuteSV --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 -mi 500 -md 500 -s 3 -t ${thread} --genotype -Ivcf ${input_vcf} -S ${sample} -L 150000 ${bam} ${reference} ${output_vcf} ${workdir}
#CCS
cueSV --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.8 -mi 500 -md 500 -s 3 -t ${thread} --genotype -Ivcf ${input_vcf} -S ${sample} -L 150000 ${bam} ${reference} ${output_vcf} ${workdir}
#ONT
cuteSV --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 -mi 500 -md 500 -s 3 -t ${thread} --genotype -Ivcf ${input_vcf} -S ${sample} -L 150000 ${bam} ${reference} ${output_vcf} ${workdir}
```

(2) LRcaller

```
#CLR/CCS/ONT
LRcaller -fa ${reference} -nt ${thread} -a seqan ${bam} ${input_vcf} ${output_vcf}
#The joint model
perl -lane 'if(/#/){print}else{print join("\t",@F[0..8],$F[12])}' ${output_vcf} > LRcaller_j_model.vcf
```

(3) Sniffles

```
#CLR/ONT
sniffles -t ${thread} -m ${bam} -v ${output_vcf} --Ivcf ${input_vcf}
#CCS
sniffles --skip_parameter_estimation -t ${thread} -m ${bam} -v ${output_vcf} --Ivcf ${input_vcf}
```

(4) SVJedi

```
#CLR
python3 svjedi.py -t ${thread} -d pb -ms 3 -v ${input_vcf} -r ${reference} -i ${input_seq_data} -o ${output_vcf}
#CCS
python3 svjedi_ccs.py -t ${thread} -d asm20 -ms 3 -v ${input_vcf} -r ${reference} -i ${input_seq_data} -o ${output_vcf}
#ONT:
python3 svjedi.py -t ${thread} -d ont -ms 3 -v ${input_vcf} -r ${reference} -i ${input_seq_data} -o ${output_vcf}
```

(5) VaPoR

```
#Vcftobed:
perl -lane 'if(/SVTYPE=INS;.*?SVLEN=([-|\d]+);.*?END=(\d+)/){$num+=1;print "$F[0]\t$F[1]\t$F[1]\tSV_$num\tINS_$F[4]"}elsif(/SVTYPE=([A-Z]+);.*?SVLEN=([-|\d]+);.*?END=(\d+)/){$num+=1;print "$F[0]\t$F[1]\t$3\tSV_$num\t$1"}' ${input_vcf} > ${output_bed}
#CLR/CCS/ONT:
vapor bed --sv-input ${input_bed} --output-path ./vapor_bed/ --output-file ${output_file} --reference ${reference} --pacbio-input ${bam}
```

**5.  Evaluation of SV genotyping methods**

```
${name}.txt: File listing paths to all variant files (on separate lines). The first line is the truth SV set, and the second line is the test SV set.

#SV merging
/The_installation_path/Jasmine-1.0.1/run.sh max_dist=1000 min_support=0 threads=2 --keep_var_ids --ignore_strand --output_genotypes file_list=${name}.txt out_file=${name}.vcf

#Generate precision, recall, and F1
SVTYPE=INS, DEL, DUP, INV, or TRA
grep "SVTYPE=$SVTYPE" ${name}.vcf|./calculate_F1.sh $SVTYPE
```

**6. Trio data**

```
${name}.txt: File listing paths to all variant files (on separate lines). Ashkenazim trio (HG002, HG003, and HG004) or a Chinese trio (HG005, HG006, and HG007).

#SV merging
/The_installation_path/Jasmine-1.0.1/run.sh max_dist=1000 min_support=0 threads=2 --keep_var_ids --ignore_strand --output_genotypes file_list=${name}.txt out_file=${name}.vcf

#Mendelian concordance
header=Mendelian_vcf_header_MFC.txt
bcftools +mendelian <(cat $header <(grep -E "SUPP_VEC=111" ${name}.vcf|grep -E -v "\.\/\."|grep -E -v "0\/0.*0\/0.*0\/0"| \
perl -lane 'if(/^[XYM]/){next}else{print}'|perl -lane '{($SVTYPE)=$F[7]=~m/SVTYPE=([^;]+)/;($LEN)=$F[7]=~m/SVLEN=([^;]+)/;($END)=$F[7]=~m/;END=([^;]+)/;}if(/GT.*?(.\/.):.*?(.\/.):.*?(.\/.):.*?/){$type="SVTYPE=".$SVTYPE;$len="SVLEN=".$LEN;$end="END=".$END;$m=$3;$f=$2;$c=$1;print join("\t",@F[0..5],"PASS",join(";",$type,$len,$end),"GT",$m,$f,$c)}'))  -t Mother,Father,Child -m c
```

**7. SV size**

```
vcf=${name}.vcf

#INS
bcftools view -i 'INFO/SVTYPE = "INS" && INFO/SVLEN >=50 && INFO/SVLEN < 100' <(sed 's/STRANDS=??;//' $vcf)|./calculate_F1.sh INS
bcftools view -i 'INFO/SVTYPE = "INS" && INFO/SVLEN >=100 && INFO/SVLEN < 1000' <(sed 's/STRANDS=??;//' $vcf)|./calculate_F1.sh INS
bcftools view -i 'INFO/SVTYPE = "INS" && INFO/SVLEN >=1000 && INFO/SVLEN < 10000' <(sed 's/STRANDS=??;//' $vcf)|./calculate_F1.sh INS
bcftools view -i 'INFO/SVTYPE = "INS" && INFO/SVLEN >=10000 && INFO/SVLEN < 1000000' <(sed 's/STRANDS=??;//' $vcf)|./calculate_F1.sh INS
Replace INS with DUP or INV. Repeat the above steps.

#DEL
bcftools view -i 'INFO/SVTYPE = "DEL" && INFO/SVLEN <= -50 && INFO/SVLEN > -100' <(sed 's/STRANDS=??;//' $vcf)|./calculate_F1.sh DEL
bcftools view -i 'INFO/SVTYPE = "DEL" && INFO/SVLEN <= -100 && INFO/SVLEN > -1000' <(sed 's/STRANDS=??;//' $vcf)|./calculate_F1.sh DEL
bcftools view -i 'INFO/SVTYPE = "DEL" && INFO/SVLEN <= -1000 && INFO/SVLEN > -10000' <(sed 's/STRANDS=??;//' $vcf)|./calculate_F1.sh DEL
bcftools view -i 'INFO/SVTYPE = "DEL" && INFO/SVLEN <= -10000 && INFO/SVLEN > -1000000' <(sed 's/STRANDS=??;//' $vcf)|./calculate_F1.sh DEL
```

**8. Tandem repeat (TR)**

```
vcf=${name}.vcf
SVTYPE=INS, DEL, DUP, INV, or TRA
TR=human_hs37d5.trf.bed

#TR
grep -E "#|SVTYPE=$SVTYPE" $vcf| bedtools intersect -a - -b $TR|sort|uniq|./calculate_F1.sh $SVTYPE

#Non-TR
grep -E "#|SVTYPE=$SVTYPE" $vcf|bedtools intersect -a - -b $TR -v|sort|uniq|./calculate_F1.sh $SVTYPE
```

**9. Imprecise breakpoint**

```
grep -E "#|SVTYPE=INS|SVTYPE=DEL" HG005_pbmm2_pbsv_5SV.vcf > HG005_pbmm2_pbsv_INS_DEL.vcf
grep -E "#|SVTYPE=DUP|SVTYPE=INV|SVTYPE=BND" HG005_pbmm2_pbsv_5SV.vcf > HG005_pbmm2_pbsv_DUP_INV_TRA.vcf

#INS and DEL
vcf=HG002_Tier1.vcf, HG002_Tier2.vcf, or HG005_pbmm2_pbsv_INS_DEL.vcf
shift_size=100, 200, 500, or 1000
python shift_INS_DEL.py hs37d5.fa $vcf ${shift_size}  ${vcf%.*}_shift_${shift_size}bp.vcf

#DUP, INV, and TRA
vcf=HG005_pbmm2_pbsv_DUP_INV_TRA.vcf
shift_size=100, 200, 500, or 1000
python shift_DUP_INV_TRA.py hs37d5.fa $vcf ${shift_size} ${vcf%.*}_shift_${shift_size}bp.vcf

Evaluate SV genotyping methods based on SVs of imprecise breakpoint.
```

