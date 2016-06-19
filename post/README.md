## Post-analysis of MIGMAP output

**IMPORTANT** MIGEC built-in analysis module requires saving a binary version of output file, which can be done by running

```bash
java -Xmx8G -jar migmap.jar -R IGH -S human --write-binary out.bin sample.fastq.gz out.txt
```

Note that ``--write-binary`` option cannot be used together with ``--by-read``.

To summarize somatic hypermutations and generate clonotype trees run

```bash
java -Xmx8G -jar migmap.jar com.antigenomics.migmap.Analyze out.bin out
```

This folder contains analysis results for IGH repertoire of hypermutating Raji cell line:

- Reads stored in ``raji_R12.fastq.gz``.. actually, those are not raw reads, but error-corrected assembled consensuses, see [MIGEC](https://github.com/mikessh/migec).
- Clonotype table and binary output are stored in ``raji_R12.txt`` and ``raji_R12.bin``.

### Analysis of somatic hyprmutation profile

![Rmd SHM analysis](https://github.com/mikessh/migmap/blob/develop/post/analyze_shm.png)

The whole analysis template is stored in R markdown format, which means it can be loaded to [Rstudio](http://rmarkdown.rstudio.com/), customized and executed. Detailed explanation of all analysis steps and plots is embedded into the template. An example PDF output is shown in ``analyze_shm.pdf`` (click [here](https://github.com/mikessh/migmap/blob/develop/post/analyze_shm.pdf) to view it).

Which will generate several text files with ``out`` prefix:

- ``out.shm.txt`` a flat file with all mutations present in sample that can be processed with ``post/analyze_shm.Rmd`` [R markdown](http://rmarkdown.rstudio.com/) template.
- ``out.net.txt``, ``out.node.txt`` and ``out.edge.txt`` (network, node and edge properties) that can be imported to [Cytoscape](http://www.cytoscape.org/) using ``Import>Table>`` menu.

For more details and an example analysis of hypermutating Raji cell repertoire go to ``post/`` folder in this repository (or click [here](https://github.com/mikessh/migmap/tree/develop/post)).