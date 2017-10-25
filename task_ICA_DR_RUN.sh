#!/bin/sh

logdir='logfiles/DR'
rm -r ${logdir}
mkdir ${logdir}
mkdir Results/DRtask_groupmaps
mkdir Results/DRtask_subjectmaps
t=1

for s in `ls -d /vols/Scratch/janineb/PROFUMO/PFM50_MSMall900/Model7_S.pfm/Subjects/??????` ; do
    S=`basename $s`
    echo $S
    echo "subID='$S'; task_ICA_DR"  >> ${logdir}/all_${t}.m
    echo "/opt/fmrib/MATLAB/R2014a/bin/matlab -singleCompThread -nojvm -nodisplay -nosplash \< ${logdir}/all_${t}.m" >> ${logdir}/AllSubmit_DR.txt
    t=`echo $t 1 + p | dc -`
done

fsl_sub -q long.q -N DRtask -t ${logdir}/AllSubmit_DR.txt -l ${logdir}