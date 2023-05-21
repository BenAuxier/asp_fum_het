
rm recA.hits.txt
echo "" >  recA.hits.txt
for i in ../../genomes/*.faa
do echo $i
blastp -query recA.faa -subject $i -outfmt 6 -evalue 1e-75 -out tmp2
cat tmp2
awk -v i=$i 'BEGIN {FS = "\t"};{OFS = "\t"}; {print i,$0}' tmp2 | sed "s/\.\.\/\.\.\/genomes\///g" >> recA.hits.txt
sleep 2s
rm tmp2
done

echo "" > recA.hits.selected.faa
#now after running to get selected hits
while read genome match protein remainder; do samtools faidx ../../genomes/$genome $protein | sed "s/>/>${genome/\.faa/}\_/g" >> recA.hits.selected.faa ; done < recA.hits.selected.txt
mafft --auto recA.hits.selected.faa > recA.hits.aligned.faa
iqtree -nt AUTO -s recA.hits.aligned.faa
