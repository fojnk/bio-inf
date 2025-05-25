package main

import (
	"fmt"
	"os"

	"github.com/scipipe/scipipe"
)

func main() {
	if len(os.Args) != 3 {
		panic("Usage: go run main.go sample.fastq.gz reference.fa")
	}
	fastqPath := os.Args[1]
	refPath := os.Args[2]

	fmt.Println("FASTQ:", fastqPath)
	fmt.Println("Reference:", refPath)

	// Создаем необходимые папки заранее
	ensureDir("fastqc_out")
	ensureDir("logs")
	ensureDir("graph")

	wf := scipipe.NewWorkflow("bwa-alignment", 4)

	// Источники
	fastqSrc := wf.NewProc("fastq_source", "")
	fastqSrc.SetOut("fastq", fastqPath)

	refSrc := wf.NewProc("ref_source", "")
	refSrc.SetOut("ref", refPath)

	// -------- FastQC --------
	// Избавляемся от вложенных директорий в названии вывода FastQC,
	// берём только базовое имя файла без пути
	fastqc := wf.NewProc("fastqc", "fastqc {i:fastq} -o fastqc_out > {o:log} 2>&1")
	fastqc.SetOut("html", "fastqc_out/{i:fastq|%.fastq.gz|basename}_fastqc.html")
	fastqc.SetOut("log", "logs/fastqc.log")
	fastqc.In("fastq").From(fastqSrc.Out("fastq"))

	// bwa index
	bwaIndex := wf.NewProc("bwa_index", "bwa index {i:ref} > {o:log} 2>&1 && touch {o:done}")
	bwaIndex.SetOut("log", "logs/bwa_index.log")
	bwaIndex.SetOut("done", "bwa_index.done")
	bwaIndex.In("ref").From(refSrc.Out("ref"))

	// bwa mem (зависит от индекса)
	bwaMem := wf.NewProc("bwa_mem", "bwa mem {i:ref} {i:reads} > {o:sam} 2> {o:log} && cat {i:index_done} > /dev/null")
	bwaMem.SetOut("sam", "aligned.sam")
	bwaMem.SetOut("log", "logs/bwa_mem.log")
	bwaMem.In("ref").From(refSrc.Out("ref"))
	bwaMem.In("reads").From(fastqSrc.Out("fastq"))

	// Создаем вход 'index_done'
	bwaMem.In("index_done").From(bwaIndex.Out("done"))

	// -------- Convert SAM to BAM --------
	sam2bam := wf.NewProc("sam2bam", "samtools view -b {i:sam} > {o:bam} 2> {o:log}")
	sam2bam.SetOut("bam", "aligned.bam")
	sam2bam.SetOut("log", "logs/sam2bam.log")
	sam2bam.In("sam").From(bwaMem.Out("sam"))

	// -------- Flagstat --------
	flagstat := wf.NewProc("flagstat", "samtools flagstat {i:bam} > {o:report} 2> {o:log}")
	flagstat.SetOut("report", "alignment_report.txt")
	flagstat.SetOut("log", "logs/flagstat.log")
	flagstat.In("bam").From(sam2bam.Out("bam"))

	// -------- Evaluate --------
	eval := wf.NewProc("evaluate", `bash -c '
	mapped_percent=$(grep -oE '[0-9]+\.[0-9]+%' -m1 {i:report} | tr -d '%')

	echo "Mapped percent: $mapped_percent%"

	if (( $(echo "$mapped_percent > 90" | bc -l) )); then
    	echo "OK"
	else
    	echo "Not OK"
	fi ' > {o:log} 2>&1`)
	eval.SetOut("log", "logs/evaluate.log")
	eval.In("report").From(flagstat.Out("report"))

	// -------- Pipeline graph --------
	wf.PlotGraphPDF("graph/result.pdf")

	// -------- Run pipeline --------
	wf.Run()
}

// ensureDir проверяет, есть ли директория, и если нет — создаёт
func ensureDir(dirName string) {
	err := os.MkdirAll(dirName, 0755)
	if err != nil {
		fmt.Printf("Error creating directory %s: %v\n", dirName, err)
		os.Exit(1)
	}
}
