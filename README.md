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
  * [ddRAD](#ddRAD)
  * [Whole-genome resequencing](#wgs)
- [Population genomic analysis](#population-genomic-analysis)
  * [Structure](#structure)
  * [DAPC](#dapc)
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

Now upload your reads.histo to [GenomeScope](http://qb.cshl.edu/genomescope/).

# Assemble with Supernova v.2

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

* BUSCOs

One metric by which to determine the completeness of your genome is by using a Benchmark of Universal Single-Copy Orthologs. These are a set of genes that are specific to certain taxa 

* Assembly statistics

** BBTools

[Download](https://sourceforge.net/projects/bbmap/files/latest/download) and [install](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/) according to the instructions.

```
export PATH=/home/ssim/SOFTWARE/BBTools/bbmap:$PATH

stats.sh in=sample1_assembly.fasta.gz >sample1_assembly.stats
```



Identifying SNPs
================

What kind of population genomic data do you have? [ddRAD](#ddRAD) or [Whole-genome resequencing](#wgs)

# ddRAD

# Whole-genome resequencing

Population genomic analysis
===========================

What kind of population genomic analysis would you like to do? [Structure](#structure) or [Discriminant analysis of principle components](#dapc)

# Structure

# Discriminant analysis of principle components

Data visualization
==================

