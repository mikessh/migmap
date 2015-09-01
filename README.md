#  HIgBlast aka IgBlastWrapper

A wrapper for [IgBlast](http://www.ncbi.nlm.nih.gov/igblast/igblast.cgi) V-(D)-J mapping tool designed to facilitate analysis immune receptor libraries profiled using high-throughput sequencing.

The software is distributed as an executable JAR file and a data bundle.

**NOTE** Last IgBlastWrp version is available [here](https://github.com/mikessh/higblast/releases/tag/v0.6) (source and readme are available [here](https://github.com/mikessh/higblast/tree/v0.6)), this is a completely re-written version of original software.

## Motivation

While being a gold standard of V-(D)-J mapping, the following limitations apply to IgBlast:

- It doesn't extract sequence of CDR3 region directly, neither provide coordinates for CDR3 region in reads. It reports reference Cys residue of Variable segment and Variable segment end in CDR3, but not Phe/Trp residue of J segment that marks the end of CDR3

- Output is not straightforward to parse and summarize, which is important to count clonotype diversity of high-throughput sequencing sample

- It doesn't account for sequence quality

- It cannot group reads into clonotypes


## Features

Present wrapper adds the following capabilities to IgBlast:

- Run on FASTQ data

- Use a comprehensive V/D/J segment database for human, mouse, rat, rabbit and rhesus monkey

- Speed-up by piping reads to IgBlast and parsing the output in parallel as the built-in ``--num-threads`` argument doesn't offer much optimization

- Assemble clonotypes, apply various filtering options such as quality filtering for CDR3 N-regions and mutations, non-coding sequence filtering, etc

- Reporting mutations (including indels) in V, D and J segments, grouped by CDR/FW region

## Pre-requisites

[Java v1.8](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html) is required to run HIgBlast. Users should then install [IgBlast v1.4.0](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone) binaries that are appropriate for their system and make sure that ``igblastn`` and ``makeblastdb`` are added to ``$PATH`` or the directory that contains binaries is specified using ``--blast-dir /path/to/bin/`` argument during HIgBlast execution.

## Execution

To see the full list of HIgBlast options run 

```bash
java -jar higblast.jar -h
```

The following command will process ``sample.fastq.gz`` file, assemble clonotypes and store them in ``out.txt``:

```bash
java -Xmx8G -jar higblast.jar -R IGH -S human sample.fastq.gz out.txt
```

HIgBlast can be also run in per-read mode and allows piping results, e.g.:

```bash
java -Xmx8G -jar higblast.jar --by-read -R IGH -S human sample.fastq.gz - | grep "IGHV1-8" > out.txt
```

## Output format

The output is provided in a tab-delimited format. Note that no header column is present in piped output. Mutations are grouped by their FW/CDR region in several columns, mutations in the same region are separated with commas. Mutation entries are stored as follows,

```
$type$position:$reference>$query
```

where ``$type`` is ``S`` for substitution, ``D`` for deletion or ``I`` for indel. Position field ``$position`` marks either the substituted base, the first deleted base or the first base after insertion. Mutation positions are provided in Variable segment coordinates with the first Variable segment germline nucleotide having position of ``0``. Reference and query bases are provided for substitution, deleted and inserted bases are provided for deletions and insertions (omitting ``>``).

Output format for assembled clonotypes is the following:

Column           | Definition
-----------------|------------------------------------------------------------------------
freq             | clonotype frequency (0..1)
count            | number of reads
mig_freq         | share of MIGs
v                | Variable segment (top hit only)
d                | Diversity segment (top hit only)
j                | Joining segment (top hit only)
cdr3nt           | CDR3 nucleotide sequence
cdr3aa           | CDR3 amino acid sequence
mutations.FR1    | Mutations in FR1 region
mutations.CDR1   | Mutations in CDR1 region
mutations.FR2    | Mutations in FR2 region
mutations.CDR2   | Mutations in CDR2 region
mutations.FR3    | Mutations in FR3 region
mutations.CDR3   | Mutations in V/D/J germline sequence in CDR3 region
mutations.FR4    | 
cdr.insert.qual  | quality string N-nucleotides in CDR3 region
mutations.qual   | mutation quality string
v.end.in.cdr3    | V segment end (exclusive) in CDR3 region
d.start.in.cdr3  | D segment start in CDR3 region or -1 if D segment not defined
d.end.in.cdr3    | D segment end (exclusive) in CDR3 region or -1 if D segment not defined
j.start.in.cdr3  | J segment start in CDR3 region or -1 if J segment not defined
v.del            | Number of nucleotides deleted from V segment 3' end
d.del.5          | Number of nucleotides deleted from D segment 5' end or -1 if D segment not defined
d.del.3          | Number of nucleotides deleted from D segment 3' end or -1 if D segment not defined
j.del            | Number of nucleotides deleted from J segment 5' end or -1 if J segment not defined
 
## Installation

See [latest release](https://github.com/mikessh/igblastwrp/releases/latest) section for HIgBlast package.

HIgBlast can be also compiled from sources using [Gradle](https://gradle.org/) with ``gradle build``. Note that in order for tests to pass IgBlast binaries should be in ``$PATH`` variable, you may need to modify following part of ``build.gradle`` 

```gradle
test {
    environment "PATH", "$System.env.PATH:/usr/local/bin/:/usr/local/ncbi/igblast/bin/"
}
```