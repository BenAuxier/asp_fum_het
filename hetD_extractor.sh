rm hetD.hits.txt
echo "" >  hetD.hits.txt

for i in ../../genomes/*.faa
do echo $i
blastp -query hetD.faa -subject $i -outfmt 6 -evalue 1e-75 -out tmp2
cat tmp2
awk -v i=$i 'BEGIN {FS = "\t"};{OFS = "\t"}; {print i,$0}' tmp2 | sed "s/\.\.\/\.\.\/genomes\///g" >> hetD.hits.txt
sleep 2s
rm tmp2
done

echo "" > hetD.hits.selected.faa
#now after running to get selected hits
while read genome match protein remainder; do samtools faidx ../../genomes/$genome $protein | sed "s/>/>${genome/\.faa/}\_/g" >> hetD.hits.selected.faa ; done < hetD.hits.selected.txt
mafft --reorder --auto hetD.hits.selected.faa > hetD.hits.aligned.faa
iqtree -nt 10 -s hetD.hits.aligned.faa
