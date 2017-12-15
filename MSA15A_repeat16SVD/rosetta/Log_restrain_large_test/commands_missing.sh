

for i in `ls -d *_0.*/`; do echo $i ; cd $i;  if [ ! -e  TPP_FL_sln.params ];then qsub -q rna rosetta.sh -N ${i:0:10} ; fi ; cd ../ ; done






