MIGMAP="java -Xmx4G -jar `ls build/libs/migmap-*.jar` --blast-dir . --data-dir data/ -S human -R IGH"
$MIGMAP post/raji_R12.fastq.gz raji.txt
$MIGMAP --by-read post/raji_R12.fastq.gz - > raji_reads.txt

for f in raji.txt raji_reads.txt
do
    if [[ ! -s $f ]]
        then exit 1
    fi
done

java -Xmx4G -cp `ls build/libs/migmap-*.jar` com.antigenomics.migmap.Analyze raji.txt

ls -lh