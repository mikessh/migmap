package igblastwrp

/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

def cli = new CliBuilder(usage: 'igblastwrp [options] input.fa[.fastq] output')
cli.h('usage')
cli.C(args: 1, argName: '\'TRA\', \'TRB\', \'TRG\', \'TRD\',  \'IGL\', \'IGK\' or \'IGH\'', 'Receptor chain [required]')
cli.S(args: 1, argName: '\'human\' or \'mouse\'', 'Species [default=HomoSapiens]')
cli.p(args: 1, 'number of threads to use [default = all available processors]')
cli.N(args: 1, 'number of reads to take')

def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() < 2 || !opt.C) {
    cli.usage()
    System.exit(0)
}

String SPECIES = opt.S ?: "human", GENE = opt.C[0..1], CHAIN = opt.C[2]
String inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1]
boolean hasD = CHAIN =~ /[BH]/

def processor = new BlastProcessor(CHAIN, hasD)

def runner = new BlastRunner(SPECIES, GENE, CHAIN, inputFileName, processor)

runner.runIdle()

//SPECIES = "human", GENE = "TR", CHAIN = "A",
//INPUT = "/Users/mikesh/Programming/igblastwrp/test_tra.fa",
//OUTPUT = "/Users/mikesh/Programming/igblastwrp/results_a.txt"
