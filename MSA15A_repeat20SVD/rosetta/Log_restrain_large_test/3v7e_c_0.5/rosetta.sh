  #!/bin/bash
  #PBS -j oe
  #PBS -l walltime=48:00:00
  #PBS -l nodes=1:ppn=1
  cd "$PBS_O_WORKDIR" 

export AMBERHOME=/cluster/apps/x86_64/packages/amber14/amber14


export PATH=$AMBERHOME/bin:$PATH:/cluster/apps/x86_64/packages/cuda/bin



RUN_CUDA="$AMBERHOME/bin/pmemd.cuda -O"

if [ ! -e  TPP_FL_sln.params ];then


sleep 13;
export ROSETTA='/home/leexa/rosetta_bin_linux_2015.39.58186_bundle';
sleep 13;
export PATH=$ROSETTA/tools/rna_tools/bin/:$PATH;
sleep 13;
source $ROSETTA/tools/rna_tools/INSTALL;
sleep 13; 
#alias python=/home/leexa/installations/Python-2.6.6/python ;

sleep 13;
export ROSETTA='/home/leexa/rosetta_bin_linux_2015.39.58186_bundle';
sleep 13;
export PATH=$ROSETTA/tools/rna_tools/bin/:$PATH;
sleep 13;
source $ROSETTA/tools/rna_tools/INSTALL;
sleep 13;
#alias python=/home/leexa/installations/Python-2.6.6/python ;

sleep 13;
export ROSETTA='/home/leexa/rosetta_bin_linux_2015.39.58186_bundle';
sleep 13;
export PATH=$ROSETTA/tools/rna_tools/bin/:$PATH;
sleep 13;
source $ROSETTA/tools/rna_tools/INSTALL;
sleep 13;
#alias python=/home/leexa/installations/Python-2.6.6/python ;



date > 1_date.txt

#cp ../*.py ./;
#python 1_automate_constain_5a_rosetta_spring.py  ;
python 1_automate_rna_helix.py  > rosetta_spring_commands;
bash rosetta_spring_commands;
source README_FARFAR;

extract_lowscore_decoys.py TPP_FL_sln.out 30;

for i in `seq 1 30` ; do echo $i;
/home/leexa/rosetta_bin_linux_2015.39.58186_bundle/main/source/bin/rna_minimize.default.linuxgccrelease -s TPP_FL_sln.out.${i}.pdb -minimize_rounds 2 -skip_coord_constraints -out:file:silent TPP2.out;done
rm TPP_FL_sln.out*pdb
date >> 1_date.txt
;fi






