#rm hetC.hits.txt

#echo "" >  hetC.hits.txt
#for i in ../../genomes/*.faa
#do echo $i
#blastp -query hetC.faa -subject $i -outfmt 6 -evalue 1e-50 -out tmp2
#cat tmp2
#awk -v i=$i 'BEGIN {FS = "\t"};{OFS = "\t"}; {print i,$0}' tmp2 | sed "s/\.\.\/\.\.\/genomes\///g" >> hetC.hits.txt
#sleep 2s
#rm tmp2
#done

echo "" > hetC.hits.selected.faa
#now after running to get selected hits
while read genome match protein remainder; do samtools faidx ../../genomes/$genome $protein | sed "s/>/>${genome/\.faa/}\_/g" >> hetC.hits.selected.faa ; done < hetC.hits.selected.txt 
mafft --reorder --auto hetC.hits.selected.faa > hetC.hits.aligned.faa
iqtree -nt 10 -s hetC.hits.aligned.faa
