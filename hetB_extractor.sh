rm hetB.hits.txt
echo "" >  hetB.hits.txt

#for i in ../../genomes/*.faa
#do echo $i
#blastp -query hetB.faa -subject $i -outfmt 6 -evalue 1e-75 -out tmp3
#cat tmp3
#awk -v i=$i 'BEGIN {FS = "\t"};{OFS = "\t"}; {print i,$0}' tmp3 | sed "s/\.\.\/\.\.\/genomes\///g" >> hetB.hits.txt
#sleep 2s
#rm tmp3
#done

echo "" > hetB.hits.selected.faa
##now after running to get selected hits
while read genome match protein remainder; do samtools faidx ../../genomes/$genome $protein | sed "s/>/>${genome/\.faa/}\_/g" >> hetB.hits.selected.faa ; done < hetB.hits.selected.txt
mafft --reorder --auto hetB.hits.selected.faa > hetB.hits.aligned.faa
iqtree -nt 10 -s hetB.hits.aligned.faa
