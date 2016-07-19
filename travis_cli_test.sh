MIGMAP="java -Xmx4G -jar `ls build/libs/migmap-*.jar`"
$MIGMAP -S human -R IGH src/test/resources/sample.fastq.gz out.txt
$MIGMAP -S human -R IGH --by-read src/test/resources/sample.fastq.gz - > out2.txt

for f in out.txt out2.txt
do
    if [[ ! -s $f ]]
        then exit 1
    fi
done