mafft-ginsi --thread 6 --auto --reorder --op 2 --ep 0.2 --adjustdirection hetA_combined.fa | tr -d "_R_" > hetA_t1_aligned.fa
mafft-ginsi --thread 6 --auto --reorder --op 3 --ep 0.2 --adjustdirection hetB_combined.fa | tr -d "_R_" > hetA_t2_aligned.fa
mafft-ginsi --thread 6 --auto --reorder --op 4 --ep 0.2 --adjustdirection hetC_combined.fa | tr -d "_R_" > hetA_t3_aligned.fa
mafft-ginsi --thread 6 --auto --reorder --op 5 --ep 0.2 --adjustdirection hetD_combined.fa | tr -d "_R_" > hetA_t4_aligned.fa
mafft-ginsi --thread 6 --auto --reorder --op 6 --ep 1.2 --adjustdirection mat_combined.fa  | tr -d "_R_" > hetA_t5_aligned.fa

echo "done going to sleep now" 
sleep 1000000m

#hetA is CM000170.1, AFUA_2G17420
#hetA is 4,650,883..4,653,737
#so to get 5kb of flanking sequence is 4,645,883..4,658,737

samtools faidx ../genomes/Afum_Af293.fna CM000170.1:4645833-4646883 | sed 's/>/>Afum_Af293_/g' > hetA.Af293.fa
samtools faidx ../genomes/Afum_Af293.fna CM000170.1:4657737-4658737 | grep -v '>' >> hetA.Af293.fa

echo "" > hetA_log.txt
echo "" > hetA_combined.fa

for i in *.fna
do echo $i
echo $i >> hetA_log.txt
blastn -max_target_seqs 1 -query hetA.Af293.fa -subject $i -evalue 1e-150 -outfmt 6 > temp_blast_results.txt
cat temp_blast_results.txt >> hetA_log.txt
MATCH_START=$(cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort -n | head -n 1)
  MATCH_END=$(cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort -n | tail -n 1)
     CONTIG=$(cat temp_blast_results.txt | cut -f 2 | head -n 1)
echo "Match is " $(expr $MATCH_END - $MATCH_START) "bases long on " $CONTIG
echo "Match start" $MATCH_START ", and match end" $MATCH_END >> hetA_log.txt
echo "Match is " $(expr $MATCH_END - $MATCH_START) "bases long on " $CONTIG >> hetA_log.txt
cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort >> hetA_log.txt
echo "" >> hetA_log.txt
if [ $((MATCH_END - MATCH_START)) -ge 9000 -a $((MATCH_END - MATCH_START)) -le 17000 ]; then samtools faidx -n 10000 $i $CONTIG:$MATCH_START-$MATCH_END | sed "s/>/>$i\_/g" >> hetA_combined.fa; fi
rm temp_blast_results.txt
done

mafft-ginsi --thread 6 --auto --reorder --adjustdirection hetA_combined.fa | tr -d "_R_" > hetA_aligned.fa

#################################################
#################################################

#hetB is CM000173.1, AFUA_5G01005
#hetB is 257,000..261,500 approximately
#so to get 5kb of flanking sequence is 252,000..266,500

samtools faidx ../genomes/Afum_Af293.fna CM000173.1:251000-252000 | sed 's/>/>Afum_Af293_/g' > hetB.Af293.fa
samtools faidx ../genomes/Afum_Af293.fna CM000173.1:265500-266500 | grep -v '>' >> hetB.Af293.fa

echo "" > hetB_log.txt
echo "" > hetB_combined.fa

for i in *.fna
do echo $i
echo $i >> hetB_log.txt
blastn -max_target_seqs 1 -query hetB.Af293.fa -subject $i -evalue 1e-150 -outfmt 6 > temp_blast_results.txt
cat temp_blast_results.txt >> hetB_log.txt
MATCH_START=$(cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort -n | head -n 1)
  MATCH_END=$(cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort -n | tail -n 1)
     CONTIG=$(cat temp_blast_results.txt | cut -f 2 | head -n 1)
echo "Match is " $(expr $MATCH_END - $MATCH_START) "bases long on " $CONTIG
echo "Match start" $MATCH_START ", and match end" $MATCH_END >> hetB_log.txt
echo "Match is " $(expr $MATCH_END - $MATCH_START) "bases long on " $CONTIG >> hetB_log.txt
cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort >> hetB_log.txt
echo "" >> hetB_log.txt
if [ $((MATCH_END - MATCH_START)) -ge 9000 -a $((MATCH_END - MATCH_START)) -le 17000 ]; then samtools faidx -n 10000 $i $CONTIG:$MATCH_START-$MATCH_END | sed "s/>/>$i\_/g" >> hetB_combined.fa; fi
rm temp_blast_results.txt
done

mafft-ginsi --thread 6 --auto --reorder --adjustdirection hetB_combined.fa | tr -d "_R_" > hetB_aligned.fa

#################################################
#################################################

#hetC is CM000174.1, AFUA_6G13810
#hetC is 3,529,780..3,530,965
#so to get 5kb of flanking sequence is 3,524,780..3,535,965

samtools faidx ../genomes/Afum_Af293.fna CM000174.1:3524780-3525780 | sed 's/>/>Afum_Af293_/g' > hetC.Af293.fa
samtools faidx ../genomes/Afum_Af293.fna CM000174.1:3534965-3535965 | grep -v '>' >> hetC.Af293.fa

echo "" > hetC_log.txt
echo "" > hetC_combined.fa

for i in *.fna
do echo $i
echo $i >> hetC_log.txt
blastn -max_target_seqs 1 -query hetC.Af293.fa -subject $i -evalue 1e-150 -outfmt 6 > temp_blast_results.txt
cat temp_blast_results.txt >> hetC_log.txt
MATCH_START=$(cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort -n | head -n 1)
  MATCH_END=$(cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort -n | tail -n 1)
     CONTIG=$(cat temp_blast_results.txt | cut -f 2 | head -n 1)
echo "Match is " $(expr $MATCH_END - $MATCH_START) "bases long on " $CONTIG
echo "Match start" $MATCH_START ", and match end" $MATCH_END >> hetC_log.txt
echo "Match is " $(expr $MATCH_END - $MATCH_START) "bases long on " $CONTIG >> hetC_log.txt
cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort >> hetC_log.txt
echo "" >> hetC_log.txt
if [ $((MATCH_END - MATCH_START)) -ge 9000 -a $((MATCH_END - MATCH_START)) -le 17000 ]; then samtools faidx -n 10000 $i $CONTIG:$MATCH_START-$MATCH_END | sed "s/>/>$i\_/g" >> hetC_combined.fa; fi
rm temp_blast_results.txt
done

mafft-ginsi --thread 6 --auto --reorder --adjustdirection hetC_combined.fa | tr -d "_R_" > hetC_aligned.fa

#################################################
#################################################

#hetD is CM000176.1, AFUA_8G01500
#hetD is 384,745..385,890
#so to get 5kb of flanking sequence is 379,745..390,890

samtools faidx ../genomes/Afum_Af293.fna CM000176.1:379745-380745 | sed 's/>/>Afum_Af293_/g' > hetD.Af293.fa
samtools faidx ../genomes/Afum_Af293.fna CM000176.1:389890-390890 | grep -v '>' >> hetD.Af293.fa

echo "" > hetD_log.txt
echo "" > hetD_combined.fa

for i in *.fna
do echo $i
echo $i >> hetD_log.txt
GENOME=$i
blastn -max_target_seqs 1 -query hetD.Af293.fa -subject $i -evalue 1e-150 -outfmt 6 > temp_blast_results.txt
cat temp_blast_results.txt >> hetD_log.txt
MATCH_START=$(cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort -n | head -n 1)
  MATCH_END=$(cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort -n | tail -n 1)
     CONTIG=$(cat temp_blast_results.txt | cut -f 2 | head -n 1)
echo "Match is " $(expr $MATCH_END - $MATCH_START) "bases long on " $CONTIG
echo "Match start" $MATCH_START ", and match end" $MATCH_END >> hetD_log.txt
echo "Match is " $(expr $MATCH_END - $MATCH_START) "bases long on " $CONTIG >> hetD_log.txt
cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort >> hetD_log.txt
echo "" >> hetD_log.txt
if [ $((MATCH_END - MATCH_START)) -ge 9000 -a $((MATCH_END - MATCH_START)) -le 17000 ]; then samtools faidx -n 10000 $i $CONTIG:$MATCH_START-$MATCH_END | sed "s/>/>$i\_/g" >> hetD_combined.fa; fi
rm temp_blast_results.txt
done

mafft-ginsi --thread 6 --auto --reorder --adjustdirection hetD_combined.fa | tr -d "_R_" > hetD_aligned.fa

#################################################
#################################################

#mat is CM000171.1, 
#AFUA_3G06160 is 1,523,105..1,524,006
#AFUA_3G06170 is 1,524,532..1,525,609

#the region of CM000171.1:1523105-1525609

#this matches p21 with region of chr3-1424808 1427315
samtools faidx ../genomes/Afum_Af293.fna CM000171.1:1518105-1519105 | sed 's/>/>Afum_Af293_/g' > mat.Af293.fa
samtools faidx ../genomes/Afum_Af293.fna CM000171.1:1529609-1530609 | grep -v '>' >> mat.Af293.fa

echo "" > mat_log.txt
echo "" > mat_combined.fa

for i in *.fna
do echo $i
echo $i >> mat_log.txt
GENOME=$i
blastn -max_target_seqs 1 -query mat.Af293.fa -subject $i -evalue 1e-150 -outfmt 6 > temp_blast_results.txt
cat temp_blast_results.txt >> mat_log.txt
MATCH_START=$(cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort -n | head -n 1)
  MATCH_END=$(cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort -n | tail -n 1)
     CONTIG=$(cat temp_blast_results.txt | cut -f 2 | head -n 1)
echo "Match is " $(expr $MATCH_END - $MATCH_START) "bases long on " $CONTIG
echo "Match start" $MATCH_START ", and match end" $MATCH_END >> mat_log.txt
echo "Match is " $(expr $MATCH_END - $MATCH_START) "bases long on " $CONTIG >> mat_log.txt
cat temp_blast_results.txt | cut -f 9-10 | tr "\t" "\n" | sort >> mat_log.txt
echo "" >> mat_log.txt
if [ $((MATCH_END - MATCH_START)) -ge 9000 -a $((MATCH_END - MATCH_START)) -le 17000 ]; then samtools faidx -n 10000 $i $CONTIG:$MATCH_START-$MATCH_END | sed "s/>/>$i\_/g" >> mat_combined.fa; fi
rm temp_blast_results.txt
done

mafft-ginsi --thread 6 --auto --reorder --adjustdirection mat_combined.fa | tr -d "_R_" > mat_aligned.fa

