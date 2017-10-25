#!/bin/sh

logdir='logfiles/CCA'
rm -r ${logdir}
mkdir ${logdir}
#rm /home/fs0/janineb/scratch/HCP/CCA/Results/CCA_all_inputs_NEW.txt
t=1

#for D in `ls /home/fs0/janineb/scratch/HCP/netmat_simulations/Results/input*using_subject_netmats.mat` ; do
#for D in `ls /home/fs0/janineb/scratch/HCP/netmat_simulations/Results/input*bin*.mat` ; do
#for D in `ls /home/fs0/janineb/scratch/HCP/netmat_simulations/Results/input_*95*.mat` ; do
for D in `ls /home/fs0/janineb/scratch/HCP/netmat_simulations/Results/input_CCA_PFM*.mat` ; do

    d=`basename $D`
    echo $d
    echo "Ddir='$d'; addpath /home/fs0/janineb/scratch/HCP/netmat_simulations/; addpath /home/fs0/janineb/scratch/HCP/CCA/files/; CCA_simulations" >> ${logdir}/all_${t}.m

    echo "matlab -singleCompThread -nojvm -nodisplay -nosplash \< ${logdir}/all_${t}.m" >> ${logdir}/AllSubmit.txt


    t=`echo $t 1 + p | dc -`

done

fsl_sub -q bigmem.q -N CCAsims -t ${logdir}/AllSubmit.txt -l ${logdir}
