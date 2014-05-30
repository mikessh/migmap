#  IgBlast wrapper  

A wrapper for [IgBlast](http://www.ncbi.nlm.nih.gov/igblast/igblast.cgi) immune repertoire analysis tool to facilitate analysis of NGS immune repertoire profiling data

The software is distributed as a bundle with .jar executable, pre-compiled platform-specific IgBlast distribution and immune gene segments library


## Motivation:

- IgBlast doesn't extract sequence of CDR3 region directly, neither provide coordinates for CDR3 region in reads. It reports reference Cys residue of Variable segment and Variable segment end in CDR3, but not Phe/Trp residue of J segment that marks the end of CDR3

- IgBlast output is not straightforward to parse and summarize, which is important to count clonotype diversity of high-throughput sequencing sample


## Features:

- Run IgBlast on FASTQ data, provide quality filtering for CDRs and alleles/hypermutations

- Speed-up by splitting reads into chunks of data and parallelizing IgBlast (the conventional --num-threads argument doesn't work well)

- Speed-up by removing redundant reads

- Report V/D/J segments and CDR1/2/3 sequences

- Group clonotypes and enumerate them

- Reporting of mismatches (alleles/hypermutations) in V gene


## Execution:

```
$> java -jar IgBlastWrapper.jar [options] inputFile outputPrefix
```

Input file could be either in FASTQ or FASTA format, raw or GZipped.

### Options:

* `-R` receptor chain. **Required**. Currently supported: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGK`, `IGL`

* `-S` species. Currently supported: `human` and `mouse`

* `-f` filter clonotypes with non-functional CDR3 region (contains stop or is out-of-frame)

* `-c` filter clonotypes with incomplete CDR3 sequence

* `-q` quality threshold. Lowest quality for CDR sequence should be higher for a clonotype to pass filter. Mutations having quality lower than the threshold are also filtered

* `-l` clonotype detalizaiton level. Possible values: `0`, `1`, `2` and `0,1`, `0,1,2`, etc. At detalization level `0` clonotypes are grouped by CDR3 sequence, all mutaitons are then assembled and enumerated within clonotype. For detalization level `1` CDR1,2 and 3 sequences are used. For level `2` CDR3 sequence and all sequence mutations are used in clonotype grouping. Output will be generated for all specified levels and `outputPrefix` will be appended with `L$level.txt`

* `-N` take `N` first reads for analysis (useful for down-sampling)

* `-p` use `p` cores (uses all cores by default)


## NOTE

IgBlastWrapper only aligns to top alleles (marked by ```*01``` in IMGT nomenclature) to speed-up. Mismatches are then extracted from alignment and reported
 
## HINT

The most straightforward way to build the scripts is to create an Intellij Project (Groovy), then use "Open Module Settings"->Artifacts(+)Jar->from modules with dependencies followed by Build->Build artifacts  