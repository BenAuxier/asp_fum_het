#bwa-mem2 index Afum_Af293.fna

#for i in fastq/*1.fastq.gz
#do echo $i
#echo ""
#SAMPLE=${i/_1.fastq.gz/}
#echo $SAMPLE
#bwa-mem2 mem -t 24 -R "@RG\tID:"$SAMPLE"\tSM:"$SAMPLE Afum_Af293.fna $SAMPLE\_1.fastq.gz $SAMPLE\_2.fastq.gz | samtools view -@ 8 -bh | samtools sort -@ 12 - > ${SAMPLE/fastq/bams}.bam
#done

#ignore outgroup populatin from Barber.2021.newick phylogeny
ls bams/*.bam | grep -v "SRR13579360" | grep -v "SRR13579435" | grep -v "SRR13579316" | grep -v "SRR13579374" | grep -v "SRR13579429" | grep -v "SRR13579436" | grep -v "SRR13579397" | grep -v "SRR13579402" | grep -v "SRR10714224" | grep -v "SRR13579463" | grep -v "SRR13579444" | grep -v "SRR13579414" > bamlist.txt



#cd bams
#for i in *.bam; do samtools index $i; done
#cd ..

#old single threaded version
freebayes -r CM000169.1 -C 15 -Z -f Afum_Af293.fna -L bamlist.txt | vcffilter -f "QUAL > 1 & SAR > 0 & SAF > 0" > barber.raw.chr1.vcf
freebayes -r CM000170.1 -C 15 -Z -f Afum_Af293.fna -L bamlist.txt | vcffilter -f "QUAL > 1 & SAR > 0 & SAF > 0" > barber.raw.chr2.vcf
freebayes -r CM000171.1 -C 15 -Z -f Afum_Af293.fna -L bamlist.txt | vcffilter -f "QUAL > 1 & SAR > 0 & SAF > 0" > barber.raw.chr3.vcf
freebayes -r CM000172.1 -C 15 -Z -f Afum_Af293.fna -L bamlist.txt | vcffilter -f "QUAL > 1 & SAR > 0 & SAF > 0" > barber.raw.chr4.vcf
freebayes -r CM000173.1 -C 15 -Z -f Afum_Af293.fna -L bamlist.txt | vcffilter -f "QUAL > 1 & SAR > 0 & SAF > 0" > barber.raw.chr5.vcf
freebayes -r CM000174.1 -C 15 -Z -f Afum_Af293.fna -L bamlist.txt | vcffilter -f "QUAL > 1 & SAR > 0 & SAF > 0" > barber.raw.chr6.vcf
freebayes -r CM000175.1 -C 15 -Z -f Afum_Af293.fna -L bamlist.txt | vcffilter -f "QUAL > 1 & SAR > 0 & SAF > 0" > barber.raw.chr7.vcf
freebayes -r CM000176.1 -C 15 -Z -f Afum_Af293.fna -L bamlist.txt | vcffilter -f "QUAL > 1 & SAR > 0 & SAF > 0" > barber.raw.chr8.vcf

/mnt/LTR_userdata/auxie001/programs/vcflib/bin/vcfcat barber.raw.chr*.vcf > barber.raw.vcf

#now using parallel
#/mnt/LTR_userdata/auxie001/programs/freebayes/scripts/freebayes-parallel <(/mnt/LTR_userdata/auxie001/programs/freebayes/scripts/fasta_generate_regions.py Afum_Af293.fna 30000) 18 -f Afum_Af293.fna -L bamlist.txt | vcffilter -f "QUAL > 1 & SAR > 0 & SAF > 0" > barber.raw.vcf
vcffilter -f "QUAL > 100 & AN > 290" barber.raw.vcf |  awk '{if (gsub(/0\/1/,/0\/1/) < 5) print $0}' | /mnt/LTR_userdata/auxie001/programs/vcflib/bin/vcfbiallelic | /mnt/LTR_userdata/auxie001/programs/vcflib/bin/vcfallelicprimitives > barber.data.vcf
