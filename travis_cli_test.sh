MIGMAP="java -Xmx4G -jar `ls build/libs/migmap-*.jar` --blast-dir . --data-dir data/ -S human -R IGH"
$MIGMAP src/test/resources/sample.fastq.gz out.txt
$MIGMAP --by-read src/test/resources/sample.fastq.gz - > out2.txt

for f in out.txt out2.txt
do
    if [[ ! -s $f ]]
        then exit 1
    fi
done