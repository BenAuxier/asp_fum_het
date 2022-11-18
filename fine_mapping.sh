#this will extract the homologous regions of the two parents used in the cross, for use both in gene prediction, as well as for synteny analysis

#Afum_p20 is the internal name for AfIR974
#Afum_p21 is the internal name ofr AfIR964
#sorry for any confusion :(

rm *.fasta *.gff3 *.paf matches.tab

samtools faidx ../genomes/Afum_p20_June5.fna chr2:4666709-4687709 > p20_hetA.fasta
samtools faidx ../genomes/Afum_p20_June5.fna chr5:0183767-0203767 > p20_hetB.fasta
samtools faidx ../genomes/Afum_p20_June5.fna chr6:3607012-3627012 > p20_hetC.fasta
samtools faidx ../genomes/Afum_p20_June5.fna chr8:0386300-0406300 > p20_hetD.fasta
samtools faidx ../genomes/Afum_p20_June5.fna chr6:1560000-1650000 > p20_hetE.fasta

echo "hetA" > matches.tab
echo "p20" >> matches.tab
minimap2 ../genomes/Afum_p20_June5.fna p20_hetA.fasta >> matches.tab
echo "p21" >> matches.tab
minimap2 ../genomes/Afum_p21_June5.fna p20_hetA.fasta >> matches.tab
echo "Af293" >> matches.tab
minimap2 GCF_000002655.1_ASM265v1_genomic.fna p20_hetA.fasta >> matches.tab
echo "" >> matches.tab

echo "hetB" >> matches.tab
echo "p20" >> matches.tab
minimap2 ../genomes/Afum_p20_June5.fna p20_hetB.fasta >> matches.tab
echo "p21" >> matches.tab
minimap2 ../genomes/Afum_p21_June5.fna p20_hetB.fasta >> matches.tab
echo "Af293" >> matches.tab
minimap2 GCF_000002655.1_ASM265v1_genomic.fna p20_hetB.fasta >> matches.tab
echo "" >> matches.tab

echo "hetC" >> matches.tab
echo "p20" >> matches.tab
minimap2 ../genomes/Afum_p20_June5.fna p20_hetC.fasta >> matches.tab
echo "p21" >> matches.tab
minimap2 ../genomes/Afum_p21_June5.fna p20_hetC.fasta >> matches.tab
echo "Af293" >> matches.tab
minimap2 GCF_000002655.1_ASM265v1_genomic.fna p20_hetC.fasta >> matches.tab
echo "" >> matches.tab

echo "hetD" >> matches.tab
echo "p20" >> matches.tab
minimap2 ../genomes/Afum_p20_June5.fna p20_hetD.fasta >> matches.tab
echo "p21" >> matches.tab
minimap2 ../genomes/Afum_p21_June5.fna p20_hetD.fasta >> matches.tab
echo "Af293" >> matches.tab
minimap2 GCF_000002655.1_ASM265v1_genomic.fna p20_hetD.fasta >> matches.tab
echo "" >> matches.tab

echo "hetE" >> matches.tab
echo "p20" >> matches.tab
minimap2 ../genomes/Afum_p20_June5.fna p20_hetE.fasta >> matches.tab
echo "p21" >> matches.tab
minimap2 ../genomes/Afum_p21_June5.fna p20_hetE.fasta >> matches.tab
echo "Af293" >> matches.tab
minimap2 GCF_000002655.1_ASM265v1_genomic.fna p20_hetE.fasta >> matches.tab
echo "" >> matches.tab

#after reading matches.tab

echo "Using values from matches.tab"
sleep 0.1

samtools faidx ../genomes/Afum_p20_June5.fna chr2:4666716-4687703 | sed "s/>/>p20_/" > p20_hetA.fasta
samtools faidx ../genomes/Afum_p21_June5.fna chr2:4678114-4699698 | sed "s/>/>p21_/" > p21_hetA.fasta
samtools faidx GCF_000002655.1_ASM265v1_genomic.fna NC_007195.1:4639576-4661152 > Af293_hetA.fasta

samtools faidx ../genomes/Afum_p21_June5.fna chr5:183768-203766 | sed "s/>/>p20_/" > p21_hetB.fasta
samtools faidx ../genomes/Afum_p20_June5.fna chr5:184268-204547 | sed "s/>/>p21_/" > p20_hetB.fasta
samtools faidx GCF_000002655.1_ASM265v1_genomic.fna NC_007198.1:246831-266778 > Af293_hetB.fasta

samtools faidx ../genomes/Afum_p20_June5.fna chr6:3607014-3627010 | sed "s/>/>p20_/" > p20_hetC.fasta
samtools faidx ../genomes/Afum_p21_June5.fna chr6:3551610-3571520 | sed "s/>/>p21_/" > p21_hetC.fasta
samtools faidx GCF_000002655.1_ASM265v1_genomic.fna NC_007199.1:3518904-3538993 > Af293_hetC.fasta

samtools faidx ../genomes/Afum_p20_June5.fna chr8:386308-406294 | sed "s/>/>p20_/" > p20_hetD.fasta
samtools faidx ../genomes/Afum_p21_June5.fna chr8:345086-365032 | sed "s/>/>p21_/" > p21_hetD.fasta
samtools faidx GCF_000002655.1_ASM265v1_genomic.fna NC_007201.1:376242-395951 > Af293_hetD.fasta

samtools faidx ../genomes/Afum_p20_June5.fna chr6:1560002-1650000 | sed "s/>/>p20_/" > p20_hetE.fasta
#remove 11kb from p21 and Af293
samtools faidx ../genomes/Afum_p21_June5.fna chr6:1566692-1664181 | sed "s/>/>p21_/" > p21_hetE.fasta
samtools faidx GCF_000002655.1_ASM265v1_genomic.fna NC_007199.1:1524915-1623571 > Af293_hetE.fasta

#get gff file from ncbi:
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.gff.gz
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.fna.gz
#for i in *.gz; do gunzip $i; done
cat p20_hetA.fasta p21_hetA.fasta > hetA.both.fasta
cat p20_hetB.fasta p21_hetB.fasta > hetB.both.fasta
cat p20_hetC.fasta p21_hetC.fasta > hetC.both.fasta
cat p20_hetD.fasta p21_hetD.fasta > hetD.both.fasta
cat p20_hetE.fasta p21_hetE.fasta > hetE.both.fasta

for i in *het*.fasta
do augustus --species=aspergillus_fumigatus --gff3=on $i > ${i/fasta/gff3}
done

minimap2 -g 100 -N 50 -p 0.1 -c p20_hetA.fasta p21_hetA.fasta > hetA.paf
minimap2 -g 100 -N 50 -p 0.1 -c p20_hetB.fasta p21_hetB.fasta > hetB.paf
minimap2 -g 100 -N 50 -p 0.1 -c p20_hetC.fasta p21_hetC.fasta > hetC.paf
minimap2 -g 100 -N 50 -p 0.1 -c p20_hetD.fasta p21_hetD.fasta > hetD.paf
minimap2 -g 100 -N 50 -p 0.1 -c p20_hetE.fasta p21_hetE.fasta > hetE.paf

tar -czvf results_fine_mapping.tar.gz *.paf *.fasta *.gff3

gffread hetA.both.gff3 -g hetA.both.fasta -y hetA.both.faa
gffread hetB.both.gff3 -g hetB.both.fasta -y hetB.both.faa
gffread hetC.both.gff3 -g hetC.both.fasta -y hetC.both.faa
gffread hetD.both.gff3 -g hetD.both.fasta -y hetD.both.faa
gffread hetE.both.gff3 -g hetE.both.fasta -y hetE.both.faa
