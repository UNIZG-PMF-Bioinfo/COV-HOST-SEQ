#!/bin/bash
# ----------------QSUB Parameters-----------------
#PBS -q q2
#PBS -M dbrkic@bioinfo.hr
#PBS -m n
#PBS -N BBmap_trimmed
#PBS -l select=ncpus=16:mem=256g:ompthreads=24
#PBS -J 0-200
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
MEMORY=256g
THREADS=16

SAMTOOLS=/samtools

BBDIR=/BBmap/bbmap
BBMAP=bbmap.sh
BBMAP_PAR="pigz=t local=t"
CHRLEN="hg38.chrom.sizes"

MEMORY=256g
THREADS=16

INPUT_DIR=/Trimming

#  input reads
IN_SEQ=($INPUT_DIR/*_1.phixtrim.fq.gz)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*${INPUT_DIR}/}
BASE=${BASE%%_*}

echo $BASE

INFILE="${FILE/_1/_#}"

$BBDIR/$BBMAP -Xmx$MEMORY threads=$THREADS in1=${INFILE} out=${BASE}.bam rgid="$BASE" rgsm="$BASE" scafstats=${BASE}_scafstats.txt mhist=${BASE}_mhist.txt bhist=${BASE}_bhist.txt qhist=${BASE}_qhist.txt aqhist=${BASE}_aqhist.txt bqhist=${BASE}_bqhist.txt ihist=${BASE}_ihist.txt ehist=${BASE}_ehist.txt qahist=${BASE}_qahist.txt indelhist=${BASE}_indelhist.txt mhist=${BASE}_mhist.txt idhist=${BASE}_idhist.txt covstats=${BASE}_covstats.txt $BBMAP_PAR &> ${BASE}.log

zcat ${BASE}.sam.gz | $SAMTOOLS view -@ $THREADS -Sb - > ${BASE}.bam
$SAMTOOLS sort -@ $THREADS ${BASE}.bam -o ${BASE}_sorted.bam
$SAMTOOLS index ${BASE}_sorted.bam
[ -f "${BASE}_sorted.bam" ] && rm -f ${BASE}.bam

$SAMTOOLS rmdup ${BASE}_sorted.bam ${BASE}_sorted_dedupe.bam
$SAMTOOLS index ${BASE}_sorted_dedupe.bam

genomeCoverageBed -ibam ${BASE}_sorted.bam -bg -split -g $CHRLEN > ${BASE}.bedGraph
wigToBigWig ${BASE}.bedGraph $CHRLEN ${BASE}.bw

[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph

