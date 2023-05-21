rm hetE.hits.txt
echo "" >  hetE.hits.txt

#for i in ../../genomes/*.faa
#do echo $i
#blastp -query hetE.faa -subject $i -outfmt 6 -evalue 1e-75 -out tmp2
#cat tmp2
#awk -v i=$i 'BEGIN {FS = "\t"};{OFS = "\t"}; {print i,$0}' tmp2 | sed "s/\.\.\/\.\.\/genomes\///g" >> hetE.hits.txt
#sleep 2s
#rm tmp2
#done

echo "" > hetE.hits.selected.faa
#now after running to get selected hits
while read genome match protein remainder; do samtools faidx -n 400 ../../genomes/$genome $protein | sed "s/>/>${genome/\.faa/}\_/g" >> hetE.hits.selected.faa ; done < hetE.hits.selected.txt
mafft --reorder --auto hetE.hits.selected.faa > hetE.hits.aligned.faa
iqtree -nt 10 -s hetE.hits.aligned.faa
