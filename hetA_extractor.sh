rm hetA.hits.txt
touch hetA.hits.txt

#for i in ../../genomes/*.faa
#do echo $i
#blastp -query hetA.faa -subject $i -outfmt 6 -evalue 1e-50 -out tmp1
#cat tmp1
#awk -v i=$i 'BEGIN {FS = "\t"};{OFS = "\t"}; {print i,$0}' tmp1 | sed "s/\.\.\/\.\.\/genomes\///g" >> hetA.hits.txt
#sleep 2s
#done

rm tmp1

echo "" > hetA.hits.selected.faa
#now after running to get selected hits
while read genome match protein remainder; do samtools faidx ../../genomes/$genome $protein | sed "s/>/>${genome/\.faa/}\_/g" >> hetA.hits.selected.faa ; done < hetA.hits.selected.txt
mafft --reorder --auto hetA.hits.selected.faa > hetA.hits.aligned.faa
iqtree -nt 10 -s hetA.hits.aligned.faa
