#!/bin/bash
set -euo pipefail

# Usage: ./pipeline.sh reads_1.fastq.gz reference.fa [NUM_THREADS]
READ=$1
REF=$2
THREADS=${3:-4}  # По умолчанию 4 потока, если не указано

echo "Running FastQC on $READ..."
fastqc -t "$THREADS" "$READ"

echo "Indexing reference genome $REF..."
bwa index "$REF"

echo "Aligning reads to reference with BWA..."
bwa mem -t "$THREADS" "$REF" "$READ" > results/aligned.sam

echo "Converting SAM to BAM..."
samtools view -@ "$THREADS" -b results/aligned.sam > results/aligned.bam

echo "Sorting BAM..."
samtools sort -@ "$THREADS" results/aligned.bam -o results/aligned_sorted.bam

echo "Indexing BAM..."
samtools index results/aligned_sorted.bam

echo "Generating alignment statistics..."
samtools flagstat results/aligned_sorted.bam > results/alignment_report.txt

mapped_percent=$(grep -oE '[0-9]+\.[0-9]+%' -m1 results/alignment_report.txt | tr -d '%')

echo "Mapped percent: $mapped_percent%"

if (( $(echo "$mapped_percent > 90" | bc -l) )); then
    echo "OK"
else
    echo "Not OK"
fi
