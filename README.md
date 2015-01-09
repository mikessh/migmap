#  IgBlast wrapper  

A wrapper for [IgBlast](http://www.ncbi.nlm.nih.gov/igblast/igblast.cgi) immune repertoire analysis tool to facilitate analysis of NGS immune repertoire profiling data

The software is distributed as a bundle with .jar executable, pre-compiled platform-specific IgBlast distribution and immune gene segments library


## Motivation:

- IgBlast doesn't extract sequence of CDR3 region directly, neither provide coordinates for CDR3 region in reads. It reports reference Cys residue of Variable segment and Variable segment end in CDR3, but not Phe/Trp residue of J segment that marks the end of CDR3

- IgBlast output is not straightforward to parse and summarize, which is important to count clonotype diversity of high-throughput sequencing sample

- IgBlast doesn't account for sequence quality


## Features:

- Run IgBlast on FASTQ data, provide quality filtering for CDRs and alleles/hypermutations

- Speed-up by splitting reads into chunks of data and parallelizing IgBlast (the conventional --num-threads argument doesn't work well)

- Speed-up by removing redundant reads

- Report V/D/J segments and CDR1/2/3 sequences

- Group clonotypes and enumerate them

- Reporting of mismatches (alleles/hypermutations) in V gene


## Execution:

```bash
java -jar path/to/igblastwrp.jar [options] inputFile outputPrefix
```

For example
```bash
java -Xmx40G -jar igblastwrp.jar -cf -R IGH --all-alleles input.fastq.gz output.txt
```

Input file could be either in FASTQ or FASTA format, raw or GZipped.

### Options:

* `-R chain` Sets the receptor chain. **Required**. Currently supported: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGK`, `IGL`

* `--all-alleles` Use all alleles during alignment (this is going to be slower). Will use only major (*01) alleles if option is not set.

* `-S species` Sets the species. Currently supported: `human`, `mouse`, `rat`, `rabbit` and `rhesus_monkey`

* `-f` filter clonotypes with non-functional CDR3 region (contains stop or is out-of-frame)

* `-c` filter clonotypes with incomplete CDR3 sequence

* `--no-cdr3` report clonotypes, for which not even a portion of CDR3 is identified

* `-q x` quality threshold. Lowest quality for CDR sequences should be higher than the threshold for a clonotype to pass filter. Mutations having quality lower than the threshold are also filtered

* `-l x` clonotype detalizaiton level. Possible values: `0`, `1`, `2` and `0,1`, `0,1,2`, etc. At detalization level `0` clonotypes are grouped by CDR3 sequence, all mutations are then assembled and enumerated within clonotype. For detalization level `1` CDR1,2 and 3 sequences are used. For level `2` CDR3 sequence and all sequence mutations are used in clonotype grouping. Output will be generated for all specified levels and `outputPrefix` will be appended with `L$level.txt`. Note that out-of-frame and stop codon presence is calculated corresponding to detalization level, i.e. stop codons in FW1 don't *noStop* field in level `2` output.

* `-N x` take `x` first reads for analysis (useful for down-sampling)

* `-p x` use `x` cores (uses all cores by default)

* `-a` [MIGEC](https://github.com/mikessh/migec) compatibility mode. Assumes FASTQ headers contain a *UMI:NNNNN:READ_COUNT* entry and performs separate read and UMI (event) counting

* `--debug` Debug mode. Will run in a single thread and pring IgBlast output to stdout. Should be used as `java -jar igblastwrp.jar --debug input_file - > results.txt`

* `-h` display help message

### Output table format:

Column       | Definition
-------------|--------------------------------------------------------------------------------------------------------
reads_count  | number of reads
reads_freq   | share of reads
mig_count    | number of MIGs (read groups with distinct UMIs), only applies when using `-a` optional
mig_freq     | share of MIGs
cdr1nt       | CDR1 nucleotide sequence
cdr2nt       | CDR2 nucleotide sequence
cdr3nt       | CDR3 nucleotide sequence
cdr1aa       | CDR1 amino acid sequence
cdr2aa       | CDR2 amino acid sequence
cdr3aa       | CDR3 amino acid sequence
inFrame      | clonotype is in-frame. Applies to CDR3 (level 0), CDR1-3 (level 1) and whole sequence (level 2)
noStop       | clonotype doesn't contain stop codons. Applies to CDR3 (level 0), CDR1-3 (level 1) and whole sequence (level 2)
complete     | CDR3 region is identified completely
vSegment     | Variable segment (major allele only)
dSegment     | Diversity segment (major allele only)
jSegment     | Joining segment (major allele only)
cdr1q        | quality string for CDR1 region
cdr2q        | quality string for CDR2 region
cdr3q        | quality string for CDR3 region
mutations    | mutations string

Mutation string is a list of mutation entries separated by "|", i.e. `entry1``|``entry2``|``...`. 

Each mutation entry is encoded in the following comma-separated format:

`reads_count:reads_freq:mig_count:mig_freq``,``nt_pos:nt_from>nt_to``,``aa_pos:aa_from>aa_to``,``region`

Here is an example:

> 3505:1E0:107:1E0,67:C>G,22:T>S,FW2|4:1,1E-3:1:9,3E-3,87:A>T,29:S>C,CDR1

## NOTE

By default **IgBlastWrapper** only aligns to top alleles (marked by ```*01``` in IMGT nomenclature) to speed-up. Mismatches are then extracted from alignment and reported.
Not all receptor & chain combinations are supported for species other than human and mouse.
 
## HINT

The most straightforward way to build the scripts is to create an Intellij Project (Groovy), then use "Open Module Settings"->Artifacts(+)Jar->from modules with dependencies followed by Build->Build artifacts. Module architecture could be changed to Maven soon..  
As for now you are recommended to use platform-specific binaries from [Latest release](https://github.com/mikessh/igblastwrp/releases/latest)