import java.nio.file.Files

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

def inputDir = new File("../igblast/db"), outputPath = "../igblast/database"
new File(outputPath).mkdir()
def id2species = ["10090": "mouse", "9606": "human"]

inputDir.listFiles().each { file ->
    def nameTokens = file.name.split("[-_\\.]")
    def species = id2species[nameTokens[1]],
        gene = nameTokens[0][0..1],
        chain = nameTokens[0][2],
        segment = nameTokens[2]

    def newName = outputPath + "/${species}_${gene}_${chain}_${segment}." + (file.name.endsWith(".fasta") ? "fa" : nameTokens[4])

    Files.copy(file.toPath(), new File(newName).toPath())
}

