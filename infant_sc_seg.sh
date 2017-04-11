#! /bin/csh

## Uses freesurfer commands to produce subcortical segmentation

## Variables

## GCA registration
mri_em_register -uns 3 -mask brainmask.mgz nu.mgz $FREESURFER_HOME/average/RB_all_2016-05-10.vc700.gca transforms/talairach.lta


mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz $FREESURFER_HOME/average/RB_all_2016-05-10.vc700.gca transforms/talairach.lta norm.mgz

mri_ca_register -align-after -nobigventricles -mask brainmask.mgz -T transforms/talairach.lta norm.mgz $FREESURFER_HOME/average/RB_all_2016-05-10.vc700.gca transforms/talairach.m3z

mri_ca_label -relabel_unlikely 9 .3 -prior 0.5 -align norm.mgz transforms/talairach.m3z $FREESURFER_HOME/average/RB_all_2016-05-10.vc700.gca aseg.auto_noCCseg.mgz

mri_cc -lta <subjid>/mri/transforms/cc_up.lta -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz <subjid>