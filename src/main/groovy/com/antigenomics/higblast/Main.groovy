/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */

package com.antigenomics.higblast

import com.antigenomics.higblast.blast.BlastInstanceFactory
import com.antigenomics.higblast.genomic.SegmentDatabase
import com.antigenomics.higblast.io.*
import com.antigenomics.higblast.mapping.ReadMappingDetailsProvider
import com.antigenomics.higblast.mapping.ReadMappingFilter

def ALLOWED_CHAINS = ["TRA", "TRB", "TRG", "TRG", "IGH", "IGL", "IGK"],
    ALLOWED_SPECIES = ["human", "mouse", "rat", "rabbit", "rhesus_monkey"],
    HOME = new File(this.class.protectionDomain.codeSource.location.path).parent.replaceAll("%20", " "),
    DEFAULT_Q = "25", ERROR_LOG = "_higblast_error.log"

def cli = new CliBuilder(usage: "higblast [options] input.(fa/fastq)[.gz] (output_file/- for stdout)")

// Runtime
cli._(longOpt: "blast-dir", args: 1, argName: "path",
        "Path to folder that contains 'igblastn' and 'makeblastdb' binaries. " +
                "[default = assume they are added to \$PATH and execute them directly]")
cli._(longOpt: "data-dir", args: 1, argName: "path",
        "Path to folder that contains data bundle (internal_data/ and optional_file/ directories). " +
                "[default = \$install_dir/data/]")
cli.n(args: 1, argName: "int",
        "Number of reads to take. [default = all]")
cli.p(args: 1, argName: "int",
        "Number of cores to use. [default = all available processors]")
cli._(longOpt: "report", args: 1, argName: "file",
        "File to store HigBlast report. Will append report line if file exists.")

// Mapping
cli.R(args: 1, argName: "chain1,...",
        "Receptor gene and chain. Several chains can be specified, separated with commas. " +
                "Allowed values: $ALLOWED_CHAINS. [required]")
cli.S(args: 1, argName: "name",
        "Species. Allowed values: $ALLOWED_SPECIES. [required]")
cli._(longOpt: "all-alleles",
        "Will use all alleles during alignment (this is going to be slower). " +
                "[default = use only major (*01) alleles]")
cli._(longOpt: "use-kabat",
        "Will use KABAT nomenclature for CDR/FW partitioning. " +
                "[default = use IMGT nomenclature]")

// Output
cli._(longOpt: "allow-incomplete",
        "Report clonotypes with partial CDR3 mapping.")
cli._(longOpt: "allow-no-cdr3",
        "Report clonotypes with no CDR3 mapping.")
cli._(longOpt: "allow-noncoding",
        "Report clonotypes that have either stop codon or frameshift in their receptor sequence.")
cli._(longOpt: "allow-noncanonical",
        "Report clonotypes that have non-canonical CDR3 (do not start with C or end with F/W residues).")
cli._(longOpt: "details", args: 1, argName: "field1,field2,.../all",
        "Additional fields to provide for output, allowed values: ${ReadMappingDetailsProvider.ALLOWED_FIELDS.join(",")}.")
cli.q(args: 1, argName: "2..40",
        "Threshold for average quality of mutations and N-regions of CDR3 [default = $DEFAULT_Q]")
cli._(longOpt: "by-read",
        "Will output mapping details for each read. " +
                "[default = assemble clonotypes and output clonotype abundance table]")
cli._(longOpt: "unmapped", args: 1, argName: "fastq[.gz]",
        "Output unmapped reads in specified file.")

// Misc
cli.h("Display this help message")

// PARSE ARGUMENTS

def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() != 2 || !opt.R || !opt.S) {
    cli.usage()
    System.exit(3)
}

// I/O

def inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1],
    reportFileName = opt.'report', unmappedFileName = opt.'unmapped'

def fastaFile = ["fasta", "fa", "fasta.gz", "fa.gz"].any { inputFileName.endsWith(it) },
    stdOutput = outputFileName == "-", byRead = (boolean) opt.'by-read'

if (!new File(inputFileName).exists()) {
    Util.error("Input file $inputFileName does not exist.", 3)
}

if (stdOutput) {
    outputFileName = null
}

def ensureDirCreated = { String fileName ->
    // Ensure all subdirs are created for output
    def parentFolder = new File(fileName).absoluteFile.parentFile
    if (parentFolder) {
        parentFolder.mkdirs()
    }
}

[outputFileName, unmappedFileName, reportFileName].each {
    if (it) {
        ensureDirCreated(it)
    }
}

// Runtime

if (opt.'blast-dir') {
    ExecutionUtil.BLAST_HOME = new File((String) opt.'blast-dir').absolutePath
}

ExecutionUtil.checkBlastBinaries()

def dataDir = new File((String) (opt.'data-dir' ?: "$HOME/data/")).absolutePath

["internal_data/", "optional_file/"].each {
    if (!new File("$dataDir/$it").exists()) {
        Util.error("Internal data folder $it missing in $dataDir.", 3)
    }
}

def threads = (opt.p ?: "$Util.N_THREADS").toInteger(), limit = (opt.n ?: "-1").toInteger()

// Mapping

def species = (String) opt.S, genes = ((String) opt.R).split(",") as List<String>,
    allAlleles = (boolean) opt.'all-alleles', useKabat = (boolean) opt.'use-kabat'

if (!ALLOWED_SPECIES.contains(species)) {
    Util.error("Unknown species $species.", 3)
}

boolean tr = false, ig = false
genes.each { gene ->
    if (!ALLOWED_CHAINS.contains(gene)) {
        Util.error("Unknown gene $gene.", 3)
    }
    if (gene.toUpperCase().startsWith("TR")) {
        tr = true
    }
    if (gene.toUpperCase().startsWith("IG")) {
        ig = true
    }
}

if (tr && ig) {
    Util.error("Mixing TCR and IG genes is not allowed.", 3)
}

// Output

def allowIncomplete = (boolean) opt.'allow-incomplete',
    allowNoCdr3 = (boolean) opt.'allow-no-cdr3', allowNoncoding = (boolean) opt.'allow-noncoding',
    allowNonCanonical = (boolean) opt.'allow-noncanonical',
    qualityThreshold = Byte.parseByte((String) (opt.q ?: DEFAULT_Q)),
    details = (opt.'details' ?: "").split(",") as List<String>

details.removeAll { it.length() == 0 }

if (!details.empty && details[0].toLowerCase() == "all") {
    details = ReadMappingDetailsProvider.ALLOWED_FIELDS
}

// RUNNING THE PIPELINE
def inputPort = fastaFile ? new FastaReader(inputFileName) : new FastqReader(inputFileName)
def outputPort = stdOutput ? StdOutput.INSTANCE : new FileOutput(outputFileName)
def detailsProvider = new ReadMappingDetailsProvider(details)
outputPort = byRead ? new ReadMappingOutput(outputPort, detailsProvider) : new ClonotypeOutput(outputPort, detailsProvider)
def blastInstanceFactory = new BlastInstanceFactory(dataDir, species, genes, allAlleles, useKabat)

if (!details.empty) {
    blastInstanceFactory.annotateV()
}

try {

    def filter = new ReadMappingFilter(qualityThreshold, allowNoCdr3, allowIncomplete, allowNoncoding, allowNonCanonical,
            unmappedFileName ? new FastqWriter(unmappedFileName) : DummyInputPort.INSTANCE)

    def pipeline = new Pipeline(inputPort, blastInstanceFactory, outputPort, filter,
            limit, threads)

    Util.report("Analyzing sample $inputFileName", 2)

    pipeline.run()

    if (reportFileName) {
        def reportFile = new File(reportFileName)

        if (!reportFile.exists()) {
            reportFile.withPrintWriter { pw ->
                pw.println("input.file\toutput.file\targuments\t" + ReadMappingFilter.OUTPUT_HEADER)
            }
        }

        reportFile.withWriterAppend { pw ->
            pw.println(inputFileName + "\t" + outputFileName + "\t" + args.join(" ") + "\t" +
                    filter.toString())
        }
    }

} catch (Exception e) {

    new File("$ERROR_LOG").withPrintWriter { writer ->
        writer.println("[${new Date()} BEGIN]")
        writer.println("[CommandLine]")
        writer.println("executing higblast.jar ${args.join(" ")}")
        writer.println("[Message]")
        writer.println(e.toString())
        writer.println("[StackTrace-Short]")
        writer.println(e.stackTrace.findAll { it.toString().contains("com.antigenomics.higblast") }.join("\n"))
        writer.println("[StackTrace-Full]")
        e.printStackTrace(new PrintWriter(writer))
        writer.println("[END]")
    }

    Util.error("Internal: ${e.toString()}. See $ERROR_LOG for details.", 1)

} finally {

    SegmentDatabase.clearTemporaryFiles()

}


