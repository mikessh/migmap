  IgBlast wrapper  
==================

A wrapper for IgBlast immune repertoire analysis tool [http://www.ncbi.nlm.nih.gov/igblast/igblast.cgi]
To be distributed as a bundle with .jar executable, pre-compiled platform-specific igblast distribution and immune gene segments library

Motivation:
- IgBlast doesn't extract sequence of CDR3 region directly, neither provide coordinates for CDR3 region in reads. It reports reference Cys residue of Variable segment and Variable segment end in CDR3, but not Phe/Trp residue of J segment that marks the end of CDR3.
- IgBlast output is not straightforward to parse and summarize, which is important to count clonotype diversity of high-throughput sequencing sample

Goals:
- Run IgBlast on FASTQ data
- Speed-up by splitting reads into chunks of data and parallelizing IgBlast (the conventional --num-threads argument doesn't work well)
- Speed-up by removing redundant reads
- Report V/D/J segments and CDR1/2/3 sequences
- Group clonotypes and enumerate them

Distant future goals:
- Reporting of mismatches (alleles/hypermutations) in V gene

NOTE: IgBlastWrapper only aligns to top alleles (marked by ```*01``` in IMGT nomenclature) to speed-up. Mismatches are then extracted from alignment and reported
 
HINT: The most straightforward way to build the scripts is to create an Intellij Project (Groovy), then use "Open Module Settings"->Artifacts(+)Jar->from modules with dependencies followed by Build->Build artifacts  