#!usr/bin/bash
# author: SiegmundWANG
# summary:
# 1. chr03_align_all_summary
# 2. inside_exon_summary
# 3. overlap_exon_summary
# 4. chr03_exon_summary




# preprocessing: we should know PRFIX!
samtools faidx genome/IRGSP-1.0_genome.fasta chr03 > genome/chr03_genome.fasta
grep "^chr03" annotation/RAP_from_pan.gff > annotation/chr03_RAP.gff
grep "^chromosome03" annotation/IRGSP4_premasked_bgf.gff > annotation/chr03_IRG.gff
grep "^Chr3" annotation/MSU_rice_genome_annotation_Release7.gff3 > annotation/chr03_MSU.gff3


# map and preprocess bam
hisat2-build genome/chr03_genome.fasta index/chr03
for file in $(ls sample/*fastq.gz); do
	hisat2 -p 32 --dta -x index/chr03 -U ${file} -S "chr03_"${file:7:9}".sam" 2> "chr03_"${file:7:9}_align_summary
done
for file in $(ls chr03*sam); do
	samtools sort -@ 32 -o ${file%.*}".bam" $file
done
for one_bam in $(ls chr03*bam); do 
	samtools index ${one_bam}
done
cat chr03*summary > chr03_align_all_summary

# merge and index new bam
samtools merge chr03_all.bam chr03*bam
samtools view -b -F 4 -@ 32 chr03_all.bam > chr03_mapped.bam
sed -i 's/Chr3/chr03/' annotation/chr03_MSU.gff3
sed -i 's/chromosome03/chr03/' annotation/chr03_IRG.gff
samtools sort -@ 32 chr03_mapped.bam -o chr03_sorted_mapped.bam
samtools index chr03_sorted_mapped.bam

# located inside exon
echo "now we go to exciting things!"
bedtools multicov -f 1 -bams chr03_sorted_mapped.bam -bed annotation/chr03_RAP.gff| grep --color=auto -v gene | awk -F "\t" '{sum+=$10;END{print sum}' > chr03_z_count_RAP &
bedtools multicov -f 1 -bams chr03_sorted_mapped.bam -bed annotation/chr03_MSU.gff3| grep --color=auto exon | awk -F "\t" '{sum+=$10;}END{print sum}' > chr03_z_count_MSU &
bedtools multicov -f 1 -bams chr03_sorted_mapped.bam -bed annotation/chr03_IRG.gff| grep --color=auto CDS | awk -F "\t" '{sum+=$10;}END{print sum}' > chr03_z_count_IRG &
cat chr03_z_count* > inside_exon_summary

# overlap exon
echo "now we are going to get overlap_exon_summary!"
bedtools multicov -bams chr03_sorted_mapped.bam -bed annotation/chr03_RAP.gff| grep --color=auto -v gene | awk -F "\t" '{sum+=$10}END{print sum}' > chr03_z2_count_RAP &
bedtools multicov -bams chr03_sorted_mapped.bam -bed annotation/chr03_MSU.gff3| grep --color=auto exon | awk -F "\t" '{sum+=$10}END{print sum}' > chr03_z2_count_MSU &
bedtools multicov -bams chr03_sorted_mapped.bam -bed annotation/chr03_IRG.gff| grep --color=auto CDS | awk -F "\t" '{sum+=$10}END{print sum}' > chr03_z2_count_IRG &
cat chr03_z2* > overlap_exon_summary

# get exon
grep --color=auto -v -w gene annotation/chr03_RAP.gff > annotation/chr03_RAP_exon.gff &
grep --color=auto -w exon annotation/chr03_MSU.gff3 > annotation/chr03_MSU_exon.gff3 &
grep --color=auto -w CDS annotation/chr03_IRG.gff > annotation/chr03_IRG_exon.gff &

bedtools multicov -bams chr03_sorted_mapped.bam -bed annotation/chr03_RAP_exon.gff| awk -F "\t" -f awk_for_RAP.awk > chr03_z3_RAP &
bedtools multicov -bams chr03_sorted_mapped.bam -bed annotation/chr03_MSU_exon.gff3| awk -F "\t" -f awk_for_MSU.awk > chr03_z3_MSU &
bedtools multicov -bams chr03_sorted_mapped.bam -bed annotation/chr03_IRG_exon.gff| awk -F "\t" -f awk_for_RAP.awk > chr03_z3_IRG &

for z3 in $(ls chr03_z3*); do
	awk -f test_validate.awk ${z3} > ${z3}_collapes
done

cat chr03_*collapse > chr03_collapse_summary

for anno in $(ls annotation/* |grep exon|grep chr03); do
	bedtools coverage -a ${anno} -b chr03_sorted_mapped.bam | awk -F "\t" -f use.awk > ${anno:11:14}_summary 
done

cat chr03*exon* > chr03_exon_summary





