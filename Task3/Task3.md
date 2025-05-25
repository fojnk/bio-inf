# Построение пайплайна получения генетических вариантов

## Стэк:
BWA + SciPipe
Работаю на MacOS

## Ссылка на загруженные прочтения из NCBI SRA

Используется рид-сет: SRR33692911 (пример WSG e.coil)
Ссылка на NCBI SRA:
https://www.ncbi.nlm.nih.gov/sra/SRR33692911

``` bash
wget -O SRR33692911.fastq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq\?acc\=SRR33692911
```

Референсный геном:
``` bash
 wget https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2/ 
```

## Bash-скрипт пайплайна 

находится в папке scripts -> pipline.sh

результат работы:
``` bash
./scripts/pipline.sh data/SRR33692911.fastq.gz data/ref_ecoli.fna.gz
Running FastQC on data/SRR33692911.fastq.gz...
application/gzip
Started analysis of SRR33692911.fastq.gz
Approx 5% complete for SRR33692911.fastq.gz
Approx 10% complete for SRR33692911.fastq.gz
Approx 15% complete for SRR33692911.fastq.gz
Approx 20% complete for SRR33692911.fastq.gz
Approx 25% complete for SRR33692911.fastq.gz
Approx 30% complete for SRR33692911.fastq.gz
Approx 35% complete for SRR33692911.fastq.gz
Approx 40% complete for SRR33692911.fastq.gz
Approx 45% complete for SRR33692911.fastq.gz
Approx 50% complete for SRR33692911.fastq.gz
Approx 55% complete for SRR33692911.fastq.gz
Approx 60% complete for SRR33692911.fastq.gz
Approx 65% complete for SRR33692911.fastq.gz
Approx 70% complete for SRR33692911.fastq.gz
Approx 75% complete for SRR33692911.fastq.gz
Approx 80% complete for SRR33692911.fastq.gz
Approx 85% complete for SRR33692911.fastq.gz
Approx 90% complete for SRR33692911.fastq.gz
Approx 95% complete for SRR33692911.fastq.gz
Analysis complete for SRR33692911.fastq.gz
Indexing reference genome data/ref_ecoli.fna.gz...
[bwa_index] Pack FASTA... 0.05 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 1.05 seconds elapse.
[bwa_index] Update BWT... 0.02 sec
[bwa_index] Pack forward-only FASTA... 0.03 sec
[bwa_index] Construct SA from BWT and Occ... 0.31 sec
[main] Version: 0.7.19-r1273
[main] CMD: bwa index data/ref_ecoli.fna.gz
[main] Real time: 1.478 sec; CPU: 1.470 sec
Aligning reads to reference with BWA...
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 169346 sequences (40000193 bp)...
[M::process] read 169548 sequences (40000193 bp)...
[M::mem_process_seqs] Processed 169346 reads in 25.092 CPU sec, 6.341 real sec
[M::process] read 169416 sequences (40000346 bp)...
[M::mem_process_seqs] Processed 169548 reads in 25.868 CPU sec, 6.308 real sec
[M::process] read 169138 sequences (40000036 bp)...
[M::mem_process_seqs] Processed 169416 reads in 27.689 CPU sec, 7.349 real sec
[M::process] read 169060 sequences (40000167 bp)...
[M::mem_process_seqs] Processed 169138 reads in 26.319 CPU sec, 6.478 real sec
[M::process] read 169086 sequences (40000273 bp)...
[M::mem_process_seqs] Processed 169060 reads in 25.474 CPU sec, 6.206 real sec
[M::process] read 142864 sequences (33831167 bp)...
[M::mem_process_seqs] Processed 169086 reads in 25.230 CPU sec, 6.196 real sec
[M::mem_process_seqs] Processed 142864 reads in 20.725 CPU sec, 5.134 real sec
[main] Version: 0.7.19-r1273
[main] CMD: bwa mem -t 4 data/ref_ecoli.fna.gz data/SRR33692911.fastq.gz
[main] Real time: 44.525 sec; CPU: 176.898 sec
Converting SAM to BAM...
Sorting BAM...
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
Indexing BAM...
Generating alignment statistics...
Mapped percent: 78.12%
Not OK
```

**Сколько не искал геном на e.coli - не смог добиться процентажа выше 78**

Пробовал брать референсный геном человека, но индексирование его с помощью **bwa** занимает слишком много времени - так и не дождался.

## Bash-скрипт для разбора samtools flagstat

Является частью скрипта выше.

## Инструкция по развертыванию и установке SciPipe на macOS

### Установка SciPipe и необходимых инструментов на macOS

Установить Homebrew (если нет):
   ```bash
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```

Установить инструменты:
``` bash
brew install bwa samtools fastqc sra-tools go graphviz
```

Установить SciPipe:

``` bash
go install github.com/scipipe/scipipe@latest
export PATH=$PATH:$(go env GOPATH)/bin
```

## Hello World на SciPipe
```
package main

import . "github.com/scipipe/scipipe"

func main() {
    wf := NewWorkflow("hello_workflow", 4)

    hello := wf.NewProc("hello", "echo 'Hello, world!' > {o:out}")
    hello.SetOut("out", "hello.txt")

    wf.Run()
}
```

Компиляция и запуск
``` bash
go run hello_world.go
```

## Пайплайн оценки качества картирования (SciPipe)

находится в папке frawework.

## Визуализация пайплайна SciPipe в виде графического файла

### Использованный способ:

SciPipe предоставляет встроенную поддержку генерации графа выполнения пайплайна в формате .dot (Graphviz). Этот файл можно преобразовать в изображение (например, .png, .svg) с помощью утилиты dot из пакета Graphviz.

FastQC-репорт, оценку выравнивания, визуализацию пайплайна и логи можно найти в папке framework

### Как работать с этим репозиторием?

 1. Есть makefile, в котором лежат сценарии запуска пайплайнов и очистки
 2. Исходные данные референсного генома и выбранного находятся в папке data
 3. Результаты работы паплайнов лежат в results (можно почистить с помощью make clean)
 4. Для замены генома надо загрузить его в data и прокинуть в makefile
