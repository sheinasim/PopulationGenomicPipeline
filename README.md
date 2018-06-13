# PopulationGenomicPipeline
Commands for performing a population genomic analysis from genome assembly to population genomic data visualization.

Table of contents
=================

- [Getting started](#getting-started)
- [Genome assembly](#genome-assembly)
  * [Estimate genome size](#estimate-genome-size)
  * [Assemble with Supernova v2](#assemble-with-supernova)
- [Assembly QA](#assembly-qc)
  * [BUSCOs](#buscos)
  * [Assembly statistics](#assembly-statistics)
- [Identifying SNPs](#identifying-snps)
- [Population genomic analysis](#population-genomic-analysis)
  * [Structure](#structure)
  * [Discriminant analysis of principle components](#dapc)
- [Data visualization](#data-visualization)

Getting started
===============

There are many considerations when embarking upon a population genomic study, each with their advantages and considerations. For example, some species already have a high-quality reference genome, but others will require a *de novo* assembly. If a reference genome is required, it is necessary to consider what type of nucleic acid samples are available as some sequencing and assembly methods require a large amount of high-molecular weight DNA or fresh material for a transcriptome assembly. Once a reference assembly is obtained, it is necessary to consider the number of samples you want to include in your study and what method you want to use to genotype your individuals.

To help with this process, here's a quick guide along with some sample commands for performing a population study.

Would you like to perform a whole-genome assembly? [Yes](#genome-assembly) or [No](#identifying-snps)

Genome assembly
===============

If you have the ability to obtain high-molecular weight DNA from a single individual of your study species, a cost-effective method for genome sequencing and assembly is a 10x linked-read. This is a method that uses high-quality Illumina barcoded short reads where all reads with the same barcode were derived from the same single-molecule of high-molecular weight DNA to infer long-range information. For genomes equal to or smaller than that of the human genome (~3gb), it is recommended to sequence a chromium library on one lane of an Illumina HiSeq X. Once that is performed, you will have three fastq files, R1, R2, and I1. These files should all be in a directory.

For this exercise, the files will be in a director called `/home/ssim/PopulationGenomicPipeline/Raw_fastq/`

## Estimate genome size

Before assembling the genome, it is necessary to know the approximate size of the genome. If you know the approximate size of your genome, [skip]() this step. 

This can be achieved using GenomeScope which estimates various aspects of your genome using a kmer count file produced by a program called jellyfish. However, to estimate kmer abundance, you must first process your fastq files using Long Ranger.

Install [Long Ranger 2.2.2](https://support.10xgenomics.com/genome-exome/software/downloads/latest).

```
## Append the directory of the program you want to use to the path so that it is available outside of that directory
export PATH=/home/ssim/SOFTWARE/longranger-2.2.2:$PATH

## Run longranger basic to process your files
longranger basic --id=sample1 --fastqs=/home/ssim/PopulationGenomicPipeline/Raw_fastq/ 
```

Install [jellyfish 2.2.10](https://github.com/gmarcais/Jellyfish/releases).

```
export PATH=/home/ssim/SOFTWARE/jellyfish-2.2.10:$PATH

## Run jellyfish 
## the .fastq file you can use to will be in sample1/outs/
jellyfish count -C -m 21 -s 1000000000 -t 10 sample1/outs/barcoded.fastq -o reads.jf
jellyfish histo -t 10 reads.jf > reads.histo
```

Now upload your reads.histo to [GenomeScope](http://qb.cshl.edu/genomescope/) and use the estimated genome size to calculate the value to use for the --maxreads flag when running Supernova v.2.

## Assemble with Supernova v.2

Install Supernova v.2 by [downloading](https://support.10xgenomics.com/de-novo-assembly/software/downloads/latest) and following [these instructions](https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/installation). 

Using your raw .fastq files use the following command with the --maxreads = (estimated genome size * 56) / 150. In this case if the estimated genome size is 370.868 Mb use: --maxreads=138457386.

```
export PATH=/home/ssim/SOFTWARE/supernova-2.0.0:$PATH

## Run Supernova 
supernova run --id=sample1 --maxreads=138457386 --fastqs=/home/ssim/PopulationGenomicPipeline/Raw_fastq/

## Write a .fastq output of your assembly only including scaffolds larger than 5kb.
supernova mkoutput --style=pseudohap --asmdir=/home/ssim/PopulationGenomicPipeline/sample1/outs/assembly --outputprefix=sample1 --minsize=5000
```

Assembly QA
===========

Once you have your reference assembly you can use various tools to assess the quality of your genome based on completeness and continuity. 

## BUSCOs

One metric by which to determine the completeness of your genome is by using a Benchmark of Universal Single-Copy Orthologs. These are a set of genes that are specific to certain taxa 

## Assembly statistics

### BBTools

[Download](https://sourceforge.net/projects/bbmap/files/latest/download) and [install](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/) according to the instructions.

```
export PATH=/home/ssim/SOFTWARE/BBTools/bbmap:$PATH

stats.sh in=sample1_assembly.fasta.gz >sample1_assembly.stats
```

Your output will look like this:
```
A       C       G       T       N       IUPAC   Other   GC      GC_stdev
0.3369  0.1630  0.1631  0.3370  0.0121  0.0000  0.0000  0.3261  0.0467

Main genome scaffold total:             7139
Main genome contig total:               9573
Main genome scaffold sequence total:    370.868 MB
Main genome contig sequence total:      366.385 MB      1.209% gap
Main genome scaffold N/L50:             31/2.738 MB
Main genome contig N/L50:               429/241.944 KB
Main genome scaffold N/L90:             185/233.255 KB
Main genome contig N/L90:               1797/35.403 KB
Max scaffold length:                    14.221 MB
Max contig length:                      1.52 MB
Number of scaffolds > 50 KB:            278
% main genome in scaffolds > 50 KB:     92.63%


Minimum         Number          Number          Total           Total           Scaffold
Scaffold        of              of              Scaffold        Contig          Contig
Length          Scaffolds       Contigs         Length          Length          Coverage
--------        --------------  --------------  --------------  --------------  --------
    All                  7,139           9,573     370,867,642     366,384,822    98.79%
   1 KB                  7,139           9,573     370,867,642     366,384,822    98.79%
 2.5 KB                  3,881           6,250     365,640,328     361,161,218    98.77%
   5 KB                  1,649           3,966     357,692,199     353,223,839    98.75%
  10 KB                    661           2,923     351,212,966     346,754,576    98.73%
  25 KB                    383           2,591     347,195,368     342,785,898    98.73%
  50 KB                    278           2,453     343,516,481     339,153,631    98.73%
 100 KB                    222           2,368     339,606,171     335,312,951    98.74%
 250 KB                    183           2,264     333,520,157     329,410,057    98.77%
 500 KB                    138           2,074     316,899,194     313,173,674    98.82%
   1 MB                     88           1,757     279,840,275     276,928,625    98.96%
 2.5 MB                     33           1,122     192,756,130     191,197,750    99.19%
   5 MB                     17             778     138,149,875     137,155,755    99.28%
  10 MB                      5             313      57,411,828      57,060,088    99.39%
```

### Assembly-stats

One way you can visualize your genome assembly is by using a program called [assembly-stats](https://github.com/rjchallis/assembly-stats). To install, clone it using git:

```
git clone https://github.com/rjchallis/assembly-stats.git
```

Then use the provided perl script to generate a .json file which will contain a summary of your genome.

```
perl asm2stats.minmaxgc.pl sample1_assembly.fasta > sample1.minmaxgc.json
```

This .json file will then be used by the assembly-stats.html which you can locally host. Append the url to get the desired [output](http://localhost/assembly-stats/assembly-stats.html?path=json/&assembly=vtam_sn2_busco&view=circle&altAssembly=vtam_hic&altView=compare&altView=cumulative&altView=table).


Repeat modeling and masking
===========================

The last step in preparing an assembly for a population genomic analysis is to mask the repeats in the genome. This is necessary so that when you align your population genomic reads to the reference, reads mapping to the repeats will not be included in the analysis. Repeat modeling and masking can be accomplished using [Repeat Modeler](http://www.repeatmasker.org/RepeatModeler/) and [Repeat Masker](http://www.repeatmasker.org/RMDownload.html). Download these and follow the instructions for installation.

```
export PATH=/home/ssim/SOFTWARE/RepeatModeler:$PATH
export PATH=/home/ssim/SOFTWARE/RepeatMasker:$PATH

# Run Repeat Modeler 
BuildDatabase -name "sample1_RM" -dir /home/ssim/PopulationGenomicPipeline/ -engine ncbi 
# This will generate a directory with a unique_name that contains a file called consensi.fa.classified

# Run Repeat Masker
RepeatMasker -pa 32 -lib unique_name/consensi.fa.classified -dir sample1_assembly_repeats_masked sample1_assembly.fasta.gz
ProcessRepeats --species drosophila sample1_assembly_repeats_masked/sample1_assembly.fasta.cat.gz
```

This will result in a file called `sample1_assembly.fasta.masked` in the `sample1_assembly_repeats_masked` directory.

Identifying SNPs
================

Now that you have a high-quality and repeat masked reference genome, you have something to align your population genomic sequences to. If you made ddRAD libraries, you can use a program called [Stacks](http://catchenlab.life.illinois.edu/stacks/) to generate many types of genotype files that can serve as input files to various analysis programs. After sequencing, your sequences will be in different files separated by index but they will require further demultiplexing by barcode. This can be achieved using the Stacks *process_radtags* and a file containing your barcode and sample information. For this exercise, your raw ddRAD sequences will be in a directory called `Raw_ddRAD`, your barcode file is `sample1_barcodes.txt`, and your demultiplexed sequences will be written to `processed_SE`.

```
# gcc is a dependency of Stacks and thus must also be loaded into the path
export PATH=/home/ssim/SOFTWARE/gcc-8_1_0/usr/local/bin:$PATH
export PATH=/home/ssim/SOFTWARE/stacks-2.0b:$PATH

# Run process_radtags to demultiplex your single-end sequences by barcode. This will result in one .fastq.gz file for each individual in the library
process_radtags -f Raw_ddRAD/sample1_R1.fastq.gz -i gzfastq -o processed_SE/ -b sample1_barcodes.txt -e nlaIII -r -c -q
```


Population genomic analysis
===========================

What kind of population genomic analysis would you like to do? [Structure](#structure) or [Discriminant analysis of principle components](#dapc)

## Structure

## Discriminant analysis of principle components

Data visualization
==================

