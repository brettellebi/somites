BAMCOUNT=/nfs/software/birney/bam-readcount/build/bin/bam-readcount
BAMPATH="../sorted_bam"
DPPATH="../dp4"
REF=../../ref/Oryzias_latipes.ASM223467v1.dna.toplevel.fa
SITES_FILE=sites.HdrR_Ho5.pos
SAMDIR=/nfs/software/birney/samtools-1.3/

FILES=($( ls $BAMPATH | grep bam | grep -v tmp | grep -v bai ))
f=${FILES[$LSB_JOBINDEX-1]}
fbname=$(basename "$f" .bam)

$SAMDIR/samtools index $BAMPATH/$f 

$BAMCOUNT -l $SITES_FILE -f $REF $BAMPATH/$f | cut -f 1,15,28,41,54,67 -d ":" | sed s/=//g | sed s/\\t:/\\t/g | sed s/:/\\t/g > $DPPATH/$fbname.dp4.txt 
