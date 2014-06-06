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

def cli = new CliBuilder(usage: 'BuildNetwork [options] input_level0 input_level1 input_level2 output_dir/')
cli.h('display help message')
def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() < 4) {
    cli.usage()
    System.exit(-1)
}

def nodeDataMap = new HashMap<String, String>(), edgeDataMap = new HashMap<String, String>()

String inputFileName1 = opt.arguments()[0], inputFileName2 = opt.arguments()[1], inputFileName3 = opt.arguments()[2],
       outputDir = opt.arguments()[3]

[inputFileName1, inputFileName2, inputFileName3].eachWithIndex { it, ind ->
    if (!new File(it).exists()) {
        println "[ERROR] Corresponding file ($it) for input level $ind does not exist."
    }
}

def HEADER = "count\tv_segment\td_segment\tj_segment\t" +
        "cdr1nt\tcdr2nt\tcdr3nt\t" +
        "cdr1aa\tcdr2aa\tcdr3aa\t" +
        "inFrame\tnoStop\tcomplete\t" +
        "cdr1q\tcdr2q\tcdr3q\t" +
        "mutations"

// todo: CDR3 hypermutations

new File(inputFileName1).splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        def nodeKey = splitLine[1..4].join("\t")
        nodeDataMap.put(nodeKey, splitLine[0..3].join("\t") + "\t" + (".\t" * 4) +
                splitLine[4..8].join("\t") + "\t" + (".\t" * 2) + splitLine[9..10].join("\t"))
    }
}

new File(inputFileName2).splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        def nodeKey = splitLine[1..6].join("\t"), upperLevelKey = splitLine[1..4].join("\t")
        nodeDataMap.put(nodeKey, splitLine.join("\t"))
        edgeDataMap.put("$nodeKey (pp) $upperLevelKey", "L1")
    }
}

new File(inputFileName3).splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        def nodeKey = splitLine[1..6].join("\t"), upperLevelKey = splitLine[[(4..6), -1].flatten()].join("\t")
        nodeDataMap.put(nodeKey, splitLine.join("\t"))
        edgeDataMap.put("$nodeKey (pp) $upperLevelKey", "L2")
    }
}

// todo: low level hypermutaions
// todo: clean by degree bottom-top

new File(outputDir).mkdirs()

//count
final static List<String> KEY_HEADER = [
        "v_segment\td_segment\tj_segment\t" +
                "cdr3nt\t" +
                "cdr3aa\t" +
                "inFrame\tnoStop\tcomplete",
        "v_segment\td_segment\tj_segment\t" +
                "cdr1nt\tcdr2nt\tcdr3nt\t" +
                "cdr1aa\tcdr2aa\tcdr3aa\t" +
                "inFrame\tnoStop\tcomplete",
        "v_segment\td_segment\tj_segment\t" +
                "cdr1nt\tcdr2nt\tcdr3nt\t" +
                "cdr1aa\tcdr2aa\tcdr3aa\t" +
                "inFrame\tnoStop\tcomplete\tmutations"
]


final static List<String> HEADE1R = ["cdr3q\tmutations",
                                     "cdr1q\tcdr2q\tcdr3q\tmutations",
                                     "cdr1q\tcdr2q\tcdr3q\tmutations"]