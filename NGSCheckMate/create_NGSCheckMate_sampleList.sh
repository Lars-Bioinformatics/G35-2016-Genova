# out=../NGSCheckMate_fastq_sampleList_T2.txt
out=$1
# echo -e "FASTQ_FILE1\tFASTQ_FILE2\tSAMPLE_NAME" > $out
for f in *R1.fastq.gz; do 
    echo -e $f"\t"${f/_R1/_R2}"\t"${f%_thruplex*} >> $out
done
