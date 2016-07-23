MIGMAP="java -Xmx4G -jar `ls build/libs/migmap-*.jar` --blast-dir . --data-dir data/ -S human -R IGH"
$MIGMAP post/raji_R12.fastq.gz raji.txt
$MIGMAP --by-read post/raji_R12.fastq.gz - > raji_reads.txt

for f in raji.txt raji_reads.txt
do
    if [[ ! -s $f ]]
        then exit 1
    fi
done

MIGMAP_MODULE="java -Xmx4G -cp `ls build/libs/migmap-*.jar`"

$MIGMAP_MODULE com.antigenomics.migmap.MergeContigs -S human -R IGH raji.txt raji_mc.txt
$MIGMAP_MODULE com.antigenomics.migmap.Correct raji_mc.txt raji_mc_corr.txt
$MIGMAP_MODULE com.antigenomics.migmap.Analyze -S human -R IGH raji_mc_corr.txt raji

ls -lh