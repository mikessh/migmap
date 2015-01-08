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

def inputFileName = "./formatted_vdjtools_refs.txt", outputPath = "./database/"
new File(outputPath).deleteDir()
new File(outputPath).mkdir()
new File("jref.txt").delete()

def speciesGeneHash = new HashSet<String>()

new File(inputFileName).splitEachLine("\t") {
    def (species, geneFull, segment, segmentFull, refPoint, seq) = it

    segmentFull = segmentFull.replaceAll("/", "_")

    speciesGeneHash.add("$outputPath/${species}_${geneFull[0..1]}_${geneFull[2]}")

    new File("$outputPath/${species}_${geneFull[0..1]}_${geneFull[2]}_${segment}.fa").withWriterAppend { writer ->
        writer.println(">$segmentFull\n$seq")
    }

    if (segment == "J") {
        new File("jref.txt").withWriterAppend { writer ->
            writer.println([species, geneFull, segmentFull, refPoint, seq].join("\t"))
        }
    }
}

speciesGeneHash.each {
    def fileName = it + "_D.fa"
    def file = new File(fileName)
    if (!file.exists()) {
        file.withPrintWriter { pw ->
            pw.println(">.\nAAAAAAAAAAAAAAAA")
        }
    }
}

new File(outputPath).listFiles().each {
    "makeblastdb -parse_seqids -dbtype nucl -in $outputPath/$it.name -out $outputPath/${it.name[0..-4]}".execute()
}