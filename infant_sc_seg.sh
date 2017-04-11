#! /bin/csh

## Uses freesurfer commands to produce subcortical segmentation in the infants

## Variables

## GCA registration: requires brainmask.mgz and nu.mgz
mri_em_register -uns 3 -mask brainmask.mgz nu.mgz $FREESURFER_HOME/average/RB_all_2016-05-10.vc700.gca transforms/talairach.lta

# CA normalization: requires ctrl_points.mgz and transforms/talairach.lta
mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz $FREESURFER_HOME/average/RB_all_2016-05-10.vc700.gca transforms/talairach.lta norm.mgz

mri_ca_register -align-after -nobigventricles -mask brainmask.mgz -T transforms/talairach.lta norm.mgz $FREESURFER_HOME/average/RB_all_2016-05-10.vc700.gca transforms/talairach.m3z

mri_ca_label -relabel_unlikely 9 .3 -prior 0.5 -align norm.mgz transforms/talairach.m3z $FREESURFER_HOME/average/RB_all_2016-05-10.vc700.gca aseg.auto_noCCseg.mgz

mri_cc -lta <subjid>/mri/transforms/cc_up.lta -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz <subjid>\



#
# recon-all --copied here for reference for debugging the above commands since I
# assume I'm going to have lots of bugs...
#
# script to run all of the reconstruction routines
#
# Original Author: Doug Greve
# CVS Revision Info:
#    $Author: zkaufman $
#    $Date: 2017/01/18 14:11:24 $
#    $Revision: 1.580.2.16 $
#
# Copyright © 2011-2016 The General Hospital Corporation (Boston, MA) "MGH"
#
# Terms and conditions for use, reproduction, distribution and contribution
# are found in the 'FreeSurfer Software License Agreement' contained
# in the file 'LICENSE' found in the FreeSurfer distribution, and here:
#
# https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
#
# Reporting: freesurfer@nmr.mgh.harvard.edu
#
#

#
umask 002;

set VERSION = '$Id: recon-all,v 1.580.2.16 2017/01/18 14:11:24 zkaufman Exp $';
set ProgName = `basename $0`;
set inputargs = ($argv);
set subjid = ();
set hemilist = (lh rh);
set mailuser = ();
set WaitForFile = ();
set nWaitsMax  = 1000;
set WaitSleep =  1m;
set NotifyFile = ();
set PrintHelp = 0;
set LF = ();
set LF_DEFAULT_NAME = recon-all.log;
set SF = ();
set SF_DEFAULT_NAME = recon-all-status.log;
set CF = ();
set CF_DEFAULT_NAME = recon-all.cmd;
set AppendLog    = 1;
set AppendStatus = 1;
set DoTime = 0;
if($?SET_FS_CMD_TIMING) set DoTime = 1;
set fs_time = "";
set ErrorFile = /dev/null
set cmd = ();

set tcsh61706 = (`tcsh --version | grep "6\.17\.06"`)
if ("$tcsh61706" != "") then
    echo ""
    echo "WARNING: tcsh v6.17.06 has an exit code bug! Please update tcsh!"
    echo ""
    # workaround to force expected behavior:
    set anyerror
endif

set Force = 0;
set DoCleanCSDF  = 0;
set DoCleanCW256 = 0;
set DoCleanTal   = 0;
set DoCleanLta   = 0;
set DoCleanCP    = 0;
set DoCleanTWM   = 0;
set DoCleanSeed  = 0;
set DoCleanPFH   = 0;
set DoCleanBM    = 0;
set DoCleanBFSE  = 0;
set DoCleanASeg  = 0;
set DoCleanWM    = 0;
set DoCleanXopts = 0;
set DoCleanT2    = 0;
set DoCleanFLAIR = 0;
set DoSuperClean = 0;

set DoShowEdits  = 0;

set InputList = ();

# parameter defaults (most can be overriden with flags)
set XOptsFile = ();
set XOptsClean = 0; # Delete a pre-existing xopts file
set XOptsUse   = 1; # Use a pre-existing xopts file (default '1': always use)
set XOptsOverwrite = 0; # Overwrite a pre-existing xopts file
set TermScriptList = ();
set ControlPointsFile = ();
set TWMControlPointsFile = ();
set PonsSeedCRS = (); # SeedCRS center of pons for mri_fill #
set CCSeedCRS = ();   # SeedCRS center of corpus callosum for mri_fill #
set RHSeedCRS = ();   # SeedCRS point in the right hemi wm for mri_fill #
set LHSeedCRS = ();   # SeedCRS point in the left  hemi wm for mri_fill #

set UseMincMritotal = 0; # if 1, then use the BIC-MNI mritotal tal reg tools
set DoTalairachUseNu = 0; # if 1, use nu.mgz as input to talairach stage
set UseYa3tTalAtlas = 0; # if 1, use 3T18yoSchwartz young-adult 3T atlas
                         # for tal_avi (recommended as 3T default by Avi S.)
set CustomTalAtlas = ""; # -custom-tal-atlas <name of atlas in average dir>
set NoEMReg = 0; # Do not run mri_em_register. Instead convert tal.xfm to lta

set DoNuMakeUchar = 1; # if 1, then run mri_nu_correct.mni using talairach.xfm
                       # to run mri_make_uchar to correct histogram problem
# Neither of these are used right now
set UseMaskNuCorrect = 0; # if 1, nu_correct runs with -mask brainmask.mgz
set DilateMaskNuCorrect = (); # grow nu_correct mask by this amount

set WaterShed = 1;      # 0=no WS, 1= WS normal, 2 = WS only, 3 = WS+1st
set WSLess    = 0;      # Shrinks skull surface
set WSMore    = 0;      # Expands skull surface
set WSPctPreFlood = (); # Pre-flooding height
set WSSeedPoint = ();   # C R S
set WSAtlas     = 0;    # 0=don't use atlas, 1=use atlas (for skull strip)
set WSGcaAtlas  = 1;    # 1=use GCA atlas and registration to do skull strip
set WSUseTalXfm = 0;    # 1=use talairach.xfm instead of talairach_with_skull
                        # to do atlas alignment during mri_watershed
set WSCopy = 0;         # Simply copy input to output ignoring other opts, 
                        # for when brain is already stripped
set DoGcut = 0;         # 1=run mri_gcut after mri_watershed
set NuIterations = 2;   # Number of iterations for nu intensity correction
set ConformMin = 0;     # 1=conformed to min dimension
set ConformKeepDC = 0;  # Keep volume direction cosines when conforming
set HiRes = 0;          # 1=hires option (conformed to min dimension)
set Norm3dIters = ();   # passed as -n to *both* mri_normalize runs
set NormMaxGrad = 1;    # passed as -g to *both* mri_normalize runs
set Norm1_b = ();       # passed as -b to the *first* mri_normalize only
set Norm1_n = ();       # passed as -n to the *first* mri_normalize only
set Norm2_b = ();       # passed as -b to the *second* mri_normalize only
set Norm2_n = ();       # passed as -n to the *second* mri_normalize only
set WMSeg_wlo = ();  # from -seg-ghi and -seg-wlo, passed to mri_segment
set WMSeg_ghi = ();  #                             and mris_make_surfaces
set FixWithGA = 1;   # for topology fixer
set FixDiagOnly = 0; # for topology fixer
set RmNeckRadius = 25;  # used by mri_remove_neck
set UseCAAlign = (-align); # flag added to mri_ca_label
set UseCAAlignAfter = (-align-after); # flag added to mri_ca_register
set UseAseg  = 1 # when set to 0 (using -noaseg),then aseg.presurf.mgz not used
                 # nor is ?h.cortex.label (which originates through the aseg)
set NoAsegInorm2 = 0 # when set to 1 (using -noaseg-inorm2), then aseg.presurf
                     # is not used during the 2nd mri_normalize step
set UseNoNeg = 0 # if 1, add '-remove_negative 1' to mris_sphere,mris_register
set NoThicken = 0 # if 1, add '-thicken 0' to mri_segment
set UnCompress = 0    # if 1, add '-uncompress' to mri_ca_reg
set BigVentricles = 0 # if 1, add '-bigventricles' to mri_ca_reg.
                      # else, add '-nobigventricles'
set DoSecondPassRenorm = 0 # if 1, add -secondpassrenorm to mri_ca_register
set UseOldTopoFix = 1 # if 1, use mris_fix_topology instead of mris_topo_fixer
set UseNewTopoFix = 0 # if 1, use mris_topo_fixer instead of mris_fix_topology
set NoRandomness = 1 # if 1, seed critical binaries with identical seeds, to
                     # ensure consistency in surface creation. otherwise, the
                     # default is to seed with current time and date,
                     # resulting in slightly different surfaces each run.
                     # affects: mris_smooth, mris_sphere, mris_topology_fixer,
                     # mris_topo_fixer, mris_ca_label,
                     # mri_robust_template (1st-base-call-only)
set RngSeed = 1234   # seed for random number generator, used only when
                     # -norandomness flag is used, and can be overriden by
                     # the flag -rng-seed <seed>
set DoMultiStrip = 0 # if 1, then multiple instances of mri_watershed and
                     # mri_em_register are run in order to determine the best
                     # skull-strip
set IsMPRAGE = 1     # if 1, then -mprage added to mri_normalize/segment, turn off with -no-mprage
set IsWashuMPRAGE = 0 # if 1, then -washu_mprage added to mri_normalize/segment
set DoConformWidth256 = 0 # if 1, then conform to 256^3 during
                          # the mri_convert -conform step
set NoNormMGZ = 0; # set to 1 when -nosubcortseg or -noaseg flag is used,
                   # which causes norm.mgz not to used during inorm2 step
set NoWMSA = 0;      # if 1, then -nowmsa flag is added to mri_ca_label
set TH3Flag = 1;   #  turns on new vertex-wise volume calc for mris_anat_stats and ?h.volume
set DoQdecCache = 0; # if 1, then create smoothed fsaverage surface files
set measurelist = ( thickness area area.pial volume curv sulc \
                    white.K white.H jacobian_white w-g.pct.mgh )
                     # ^ these are the files smoothed by DoQdecCache (-qcache)
set UserMeasureList = 0; # if 1, then 'measurelist' gets -measure args
set measuredir = (); # for specifying an alternate path to measure files
set fwhmlist = ( 0 5 10 15 20 25 ) # see DoQdecCache
set target = fsaverage # see DoQdecCache
set SmoothCortexOnly = 1; # For -qcache. smooths only within ?h.cortex.label

set DoMakefile = 0; # if 1, run make -f recon-all.makefile $(MakefileTarget)
set MakefileTarget = () # argument to recon-all -make
set DoLabelV1 = 0; # if 1, create V1 label from O.Hinds V1 prediction atlas
set DoRobustMotionCor = 1; # if 1, then use mri_robust_template for motion cor
set mc_robust_template_avg_arg = 1; # when using mri_robust_template for motion
                     # correction, construct template from: 0 Mean, 1 Median
set UseCuda = 0; # if 1 (-use-cuda), then use GPU versions of tools:
                 # mri_em_register_cuda, mri_ca_register_cuda,
                 # mris_inflate_cuda, mris_sphere_cuda
set GetCuda = 0; # if 1 (-get-cuda) then print cuda info and exit
set PialNoAparc = 0; # if 1 (-pial-noaparc), then add -noaparc flag to
                     # mris_make_surfaces to bypass usage of parcellation
set UseCPsWithCaNorm = 0 # -canorm-usecps enable control points with ca_norm
set DoT2pial = 0; # if 1, mris_make_surfaces refines pial using T2
set T2max = ();
set DoFLAIRpial = 0; # if 1, mris_make_surfaces refines pial using FLAIR

# set multi-threaded stuff to single-threaded/cpu for cluster politeness.
# various binaries (em_reg, ca_reg, mris_sphere) have OpenMP enabled code.
# -openmp <num_threads> allows specifying more threads.
# to support -make, where recon-all is re-entrant, look for FS_OMP_NUM_THREADS
if ($?FS_OMP_NUM_THREADS) then
  setenv OMP_NUM_THREADS $FS_OMP_NUM_THREADS # var declared in -openmp
else
  setenv OMP_NUM_THREADS 1 # default is one thread, -openmp will override
endif
set OMP_NUM_SET = 0 # set to 1 if -openmp <num> is used

# hippocampal subfields processing uses ITK threading
setenv ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS 1 # -itkthreads <num_threads>
set DoParallel = 0; # if 1, then run with -openmp 4 and -itkthreads 4, 
                    # and run rh/lh stages in parallel

# Longitudinal processing:
set longitudinal = 0;    # if 1, will run the longitudinal scheme
set longbaseid = ();
set tpNid = ();
set longbasetotpN_regfile = ();   # reg file to align longbase to current subj
set UseConcatLtaTal = 0; # if 1, use mri_concatenate_lta during tal creation
set UseLongbaseCtrlVol = 0;   # if 1, use ctrl volume of longbase in norm step
set UseLongbaseWMedits = 0;   # if 1, xfer wm edits from base (default is cross)
set UseAsegFusion = 1; # if 0, dont create 'fused' asegs from timepoints
set DoCreateBaseSubj = 0; # if 1, create the 'base' subject used in longitud
set BaseSubjInvol = (orig.mgz); # -base-invol allows using some other file
set BaseSubjsList = (); # subject set grabbed from -base-tp args
set BaseSubjsListFname = (base-tps); # file containing BaseSubjsList
set robust_template_avg_arg = 1; # construct template from: 0 Mean, 1 Median
set DoNuIntensityCor3T = 0; # if 1, use Zheng, Chee, Zagorodnov params for 3T
set DoAddTp = 0; # if 1, then 'fake'-add this timepoint to long subj set
set DoAffineBase = 0; # if 1, allow affine when creating base (fix calibration)
set UseCubic = 0; # if 1, use cubic spline when mri_convert does conform step
set UseFixMtl = 0; # if 1, add -fix_mtl flag to make_surfaces

# Hippocampal subfields:
set DoHippoSF_T1 = 0;  # if 1, run hippo subfield seg, exvivo atlas (only T1)
set DoHippoSF_T1T2 = 0;# if 1, run hippo subfield seg, exvivo atlas (T1 and T2)
set DoHippoSF_T2 = 0;  # if 1, run hippo subfield seg, exvivo atlas (only T2)
set T2volHippoSF = ();

# Brainstem substructures
set DoBSsubst = 0;

# For defacing, as found in:
# $FREESURFER_HOME/average/
set brain_template = talairach_mixed_with_skull.gca
set face_template  = face.gca

# For subcortical segmentation
set GCA      = RB_all_2016-05-10.vc700.gca
set GCASkull = RB_all_withskull_2016-05-10.vc700.gca
set GCADIR   = "${FREESURFER_HOME}/average"
set GCAOutputName = ();
set GCARegIterations = ();
set GCARegTol  = ();

# For cortical registration, as found in $AvgCurvTifPath/$hemi.$AvgCurvTif
set AvgCurvTifPath = "${FREESURFER_HOME}/average"
#set AvgCurvTif = average.curvature.filled.buckner40.tif
#set AvgCurvTif = curvature.buckner40.2016-03-20.tif
set AvgCurvTif = folding.atlas.acfb40.noaparc.i12.2016-08-02.tif

# Desikan-Killiany cortical parcellation atlas (-cortparc), as found in:
# $FREESURFER_HOME/average/$hemi.$GCS
# The 2009-03-04 atlas contains the insula label.
# The 2010-03-25 atlas has a different color for temporalpole (the old
# color was gray, which looked like the default tksurfer surface color).
set OLD_OLD_GCS = curvature.buckner40.filled.desikan_killiany.2007-06-20.gcs
set OLD_GCS = curvature.buckner40.filled.desikan_killiany.2009-03-04.gcs
set NEW_OLD_GCS = curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs
# The GCSs must have been created with the matching folding atlas
set GCS = DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
set GCSDIR = "${FREESURFER_HOME}/average"

# Christophe Destrieux cortical parcellation atlas (-cortparc2):
set OLD_OLD_DESTRIEUX_GCS = atlas2002_simple.gcs
set OLD_OLD_DESTRIEUX_NAME = a2002s
set OLD_DESTRIEUX_GCS = atlas2005_simple.gcs
set OLD_DESTRIEUX_NAME = a2005s
set NEW_OLD_DESTRIEUX_GCS = destrieux.simple.2009-07-29.gcs
set DESTRIEUX_GCS = CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
set DESTRIEUX_NAME = a2009s

# Mindboggle cortical parcellation atlas (-cortparc3):
set DKTATLAS40_GCS = DKTatlas40.gcs
set DKTATLAS40_NAME = DKTatlas40
set DKTATLAS_GCS = DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
set DKTATLAS_NAME = DKTatlas

# Initialization method for BBR
set BBRInit = coreg; # others are spm, header, rr, best

# Set this to 0 to do everything but run the command.
# This is good for debugging.
set RunIt = 1;

set PatchDir = ();

# print versions and exit
set DoVersionsOnly = 0;

#----- Volume -----------#
set DoConvertInput   = 0;
set DoConvertT2Input = 0;
set DoConvertFlairInput = 0;
set DoCreateBaseInput = 0;
set DoMotionCor      = 0;
set DoDeface         = 0;
set DoNuIntensityCor = 0;
set DoTalairach      = 0;
set DoTalCheck       = 0;
set DoNormalization  = 0;
set DoNormalization2 = 0;
set DoMaskBFS        = 0;
set UseControlPoints = 0;
set UseTWMControlPoints = 0;
set DoSkullStrip     = 0;
set DoSegmentation   = 0;
set DoGCAReg         = 0;
set DoCARegInv       = 0;
set DoCANormalize    = 0;
set DoCAReg          = 0;
set DoRemoveNeck     = 0;
set DoSkullLTA       = 0;
set DoCALabel        = 0;
set DoASegMerge      = 0;
set DoFill           = 0;
#----- Surface -----------#
set DoTessellate     = 0;
set SvInitOrigSurf   = 0;
set DoSmooth1        = 0;
set DoInflate1       = 0;
set DoQSphere        = 0;
set DoFix            = 0;
set DoSmooth2        = 0;
set DoInflate2       = 0;
set DoCurvHK         = 0;
set DoCurvStats      = 0;
set DoSphere         = 0;
set DoSurfReg        = 0;
set SurfRegToSubj    = ();
set DoJacobianWhite  = 0;
set DoJacobianDist0  = 0;
set DoContraSurfReg  = 0;
set DoContraSurfRegWithinSubject = 0;
set DoAvgCurv        = 0;
set DoMorphRGB       = 0;
set DoWhiteSurfs     = 0;
set DoCortParc       = 0;
set DoCortParc2      = 0;
set DoCortParc3      = 0;
set DoPialSurfs      = 0;
set DoPctSurfCon     = 0;
set DoSurfVolume     = 0;
set DoParcStats      = 0;
set DoParcStats2     = 0;
set DoParcStats3     = 0;
set DoLocalGyriIndex = 0;
set DoBaLabels       = 0;
set DoLabelExvivoEC  = 0;
# ----------- Surface and Volume -----------#
set DoCortRibbonVolMask = 0;
set DoRelabelHypos   = 0;
set DoAParc2ASeg     = 0;
set DoAPas2ASeg      = 0;
set DoSegStats       = 0;
set DoWMParc         = 0;
set DoAParcASegStats = 0;

set DoVnoMatchCheck  = 0;

set DoIsRunning  = 1;
set IsRunningFile  = ();
set DoneFile = ();

setenv LANG C # Required by MNI tool

# -------------------------------------------------- #
set argv0 = ($argv); # make a copy
set PWD = pwd;
# better yet, make sure the real pwd is used:
if ( -e /bin/pwd ) set PWD = /bin/pwd
if($#argv == 0) goto usage_exit;
set n = `echo $argv | egrep -e -help | wc -l`
if($n != 0) then
  set PrintHelp = 1;
  goto usage_exit;
endif
set n = `echo $argv | egrep -e -version | wc -l`
if($n != 0) then
  cat $FREESURFER_HOME/build-stamp.txt
  exit 0;
endif

source $FREESURFER_HOME/sources.csh

goto parse_args;
parse_args_return:
goto check_params;
check_params_return:

set StartTime = `date`;
set tSecStart = `date '+%s'`;

# This allows the user to require that the build stamp be
# consistent from one recon-all invocation to the next.
# Good for frozen versions.
if($?REQUIRE_FS_MATCH == 0) setenv REQUIRE_FS_MATCH 0
#echo "REQUIRE_FS_MATCH $REQUIRE_FS_MATCH"
set bstampfile0 = $FREESURFER_HOME/build-stamp.txt
mkdir -p $SUBJECTS_DIR/$subjid/scripts
set bstampfile  = $SUBJECTS_DIR/$subjid/scripts/build-stamp.txt
if(-e $bstampfile0) then
  if(! -e $bstampfile) cp $bstampfile0 $bstampfile
  set bstamp0 = `cat $bstampfile0`
  set bstamp  = `cat $bstampfile`
  if("$bstamp0" != "$bstamp") then
    if($REQUIRE_FS_MATCH) then
      echo "ERROR: FreeSurfer build stamps do not match"
      echo "Subject Stamp: $bstamp"
      echo "Current Stamp: $bstamp0"
      exit 1;
    else
      echo "INFO: FreeSurfer build stamps do not match"
    endif
  endif
  echo "Subject Stamp: $bstamp"
  echo "Current Stamp: $bstamp0"
endif
rm -f $SUBJECTS_DIR/$subjid/scripts/patchdir.txt
if($?PatchDir) then
  echo "$PatchDir" > $SUBJECTS_DIR/$subjid/scripts/patchdir.txt
endif
cp $bstampfile0 $SUBJECTS_DIR/$subjid/scripts/lastcall.build-stamp.txt

if ($DoMakefile) then
  setenv RECONALL_MAKE_SUBJECT $subjid
  set make_flags=( $MakefileTarget )
  if ( ! $RunIt) set make_flags=( -n $make_flags )
  echo "Subject '$subjid': make $make_flags"
  make -f $FREESURFER_HOME/bin/recon-all.makefile ${make_flags}
  set makestatus=($status)
  unsetenv RECONALL_MAKE_SUBJECT
  exit ($makestatus)
endif

if ($DoTime) then
  fs_time ls >& /dev/null
  if ( ! $status) set fs_time=(fs_time)
endif

echo "INFO: SUBJECTS_DIR is $SUBJECTS_DIR"

# Get "True" FS HOME
pushd $FREESURFER_HOME > /dev/null
set freesurfer_home_true = `pwd`;
popd > /dev/null
echo "Actual FREESURFER_HOME $freesurfer_home_true"

set DateString = "`date '+%y%m%d%H%M'`"

# -superclean deletes everything except contents of mri/orig
if($DoSuperClean) then
  echo "super-cleaning..."
  set cmd = (rm -Rfv \
    $subjdir/label \
    $subjdir/scripts \
    $subjdir/stats \
    $subjdir/surf \
    $subjdir/tmp \
    $subjdir/touch \
    $subjdir/trash \
    $subjdir/mri/transforms)
  echo "$cmd"
  if($RunIt) $cmd
  set cmd = (rm -v $subjdir/mri/* $subjdir/mri/.xdebug*)
  echo "$cmd"
  if($RunIt) $cmd
endif

cd $subjdir # This variable is set in check_params
mkdir -p mri scripts surf tmp label touch stats touch trash
mkdir -p mri/transforms mri/transforms/bak mri/orig
set touchdir = $subjdir/touch

# Create cmd and env files from scratch
if(! $DoVersionsOnly) then
  set CF = ($subjdir/scripts/$CF_DEFAULT_NAME)
  # rm -f $CF, let it accumulate commands
  echo "\n\n#---------------------------------" >> $CF
  echo "# New invocation of recon-all `date` " >> $CF
  # Create a separate file for the env
  set ENVF = $subjdir/scripts/recon-all.env
  if(-e $ENVF) mv -f $ENVF $ENVF.bak
  date                                     >> $ENVF
  echo "FREESURFER_HOME $FREESURFER_HOME"  >> $ENVF
  echo "Actual FREESURFER_HOME $freesurfer_home_true"  >> $ENVF
  pwd                                      >> $ENVF
  echo "setenv SUBJECTS_DIR $SUBJECTS_DIR" >> $ENVF
  echo $inputargs                          >> $ENVF
  uname -a                                 >> $ENVF
  echo ""                                  >> $ENVF
  limit                                    >> $ENVF
  echo ""                                  >> $ENVF
  printenv                                 >> $ENVF
endif

if($DoVersionsOnly) then
  if (-e /dev/stdout) then
    set LF = /dev/stdout
    set SF = /dev/stdout
  else
    set LF = /dev/null
    set SF = /dev/null
  endif
endif

# ------------ Create the log file --------------- #
if($#LF == 0) then
  set LF = ($subjdir/scripts/$LF_DEFAULT_NAME)
  if(-e $LF) then
    ls -l $LF
    if(! $AppendLog) then
      mv -f $LF $LF.old
    else
      # if running using -make, then dont bother with repeated info dumps
      if ($?RECONALL_MAKE_SUBJECT) goto skip_new_invo
      echo "\n\n"  >> $LF
      echo "New invocation of recon-all "  >> $LF
      echo "\n\n"  >> $LF
    endif
  endif
else
  if(-e $LF) then
    if ($?RECONALL_MAKE_SUBJECT) goto skip_new_invo
    echo "\n\n"  >> $LF
    echo "New invocation of recon-all "  >> $LF
    echo "\n\n"  >> $LF
  endif
endif
skip_new_invo:

date >> $LF
$PWD >> $LF
echo $0 >> $LF
echo $inputargs >> $LF
# if running using -make, then dont bother with repeated info dumps
if ($?RECONALL_MAKE_SUBJECT) goto skip_all_info
echo "subjid $subjid" >> $LF
echo "setenv SUBJECTS_DIR $SUBJECTS_DIR" >> $LF
echo "FREESURFER_HOME $FREESURFER_HOME" >> $LF
echo "Actual FREESURFER_HOME $freesurfer_home_true" >> $LF
if (-e $FREESURFER_HOME/build-stamp.txt) then
  echo "build-stamp.txt: `cat $FREESURFER_HOME/build-stamp.txt`" >> $LF
endif
uname -a | tee -a $LF
limit >> $LF
if (-e /usr/bin/free) then
  echo "" >> $LF
  /usr/bin/free >> $LF
  echo "" >> $LF
endif
if ("`uname -s`" == "Darwin") then
  echo "" >> $LF
  /usr/bin/top -l 1 | grep PhysMem >> $LF
  echo "" >> $LF
endif
if($?PBS_JOBID) then
  # If the job has been submitted to launchpad, get jobid
  echo "pbsjob $PBS_JOBID"  >> $LF
endif

# check for existence of bc (binary calculator)                                                                                               
# some minimal installs of centos do not have it
set cmd = (which bc)
$cmd >& /dev/null
if($status) then
  echo "ERROR: OS is missing bc (binary calculator) utility" |& tee -a $LF
  exit 1;
endif

## gather all versions here
echo "########################################" >> $LF
echo "program versions used" >> $LF
echo $VERSION                >> $LF
mri_motion_correct.fsl -version >> $LF
if (-e $FREESURFER_HOME/bin/flirt.fsl) flirt.fsl -version >> $LF
talairach_avi --version >> $LF
tkregister2_cmdl --all-info >> $LF
nu_correct -version >> $LF
mri_make_uchar -all-info >> $LF
mri_normalize -all-info >> $LF
mri_watershed -all-info >> $LF
mri_gcut -all-info >> $LF
mri_segment -all-info >> $LF
mri_label2label -all-info >> $LF
mri_em_register -all-info >> $LF
mri_ca_normalize -all-info >> $LF
mri_ca_register -all-info >> $LF
mri_ca_label -all-info >> $LF
mri_pretess -all-info >> $LF
mri_fill -all-info >> $LF
mri_tessellate -all-info >> $LF
mri_concatenate_lta -all-info >> $LF
mri_normalize_tp2 -all-info >> $LF
mris_smooth -all-info >> $LF
mris_inflate -all-info >> $LF
mris_curvature -all-info >> $LF
mris_sphere -all-info >> $LF
mris_fix_topology -all-info >> $LF
mris_topo_fixer -all-info >> $LF
mris_ca_label -all-info >> $LF
mris_euler_number -all-info >> $LF
mris_make_surfaces -all-info >> $LF
mris_register -all-info >> $LF
mris_volmask --all-info >> $LF
mris_anatomical_stats -all-info >> $LF
mrisp_paint -all-info >> $LF
mris_curvature_stats -all-info >> $LF
if(-e .xdebug_mris_curvature_stats) rm -f .xdebug_mris_curvature_stats
mris_calc -all-info >> $LF
if(-e .xdebug_mris_calc) rm -f .xdebug_mris_calc
mri_robust_register -all-info >> $LF
mri_robust_template -all-info >> $LF
mri_and -all-info >> $LF
mri_or -all-info >> $LF
mri_fuse_segmentations -all-info >> $LF
mri_segstats -all-info >> $LF
mri_relabel_hypointensities -all-info >> $LF

echo "#######################################" >> $LF
echo "GCADIR $GCADIR" >> $LF
echo "GCA $GCA" >> $LF
echo "GCASkull $GCASkull" >> $LF
echo "AvgCurvTif $AvgCurvTif" >> $LF
echo "GCSDIR $GCSDIR" >> $LF
echo "GCS $GCS" >> $LF
echo "#######################################" >> $LF
skip_all_info:

if($DoVersionsOnly) exit 0;

# Delete the error file. This is created when error_exit is run.
set ErrorFile = $subjdir/scripts/recon-all.error
rm -f $ErrorFile
# Delete the done file. This is created when recon-all exits normally
if($#DoneFile == 0) then
  set DoneFile = $subjdir/scripts/recon-all.done
  rm -f $DoneFile
endif

# ------------ Create the IsRunning File --------- #
if($DoIsRunning) then
  set IsRunningLH   = $subjdir/scripts/IsRunning.lh
  set IsRunningRH   = $subjdir/scripts/IsRunning.rh
  set IsRunningLHRH = $subjdir/scripts/IsRunning.lh+rh
  set bailfile = ();
  if($#hemilist == 1) then
    set hemi = $hemilist;
    set IsRunningFile = $subjdir/scripts/IsRunning.$hemi
    if(-e $IsRunningLHRH) set bailfile = $IsRunningLHRH
  else
    set IsRunningFile = $subjdir/scripts/IsRunning.lh+rh
    if(-e $IsRunningLH)   set bailfile = $IsRunningLH
    if(-e $IsRunningRH)   set bailfile = $IsRunningRH
  endif
  if(-e $IsRunningFile) set bailfile = $IsRunningFile
  if($#bailfile) then
    echo ""
    echo "ERROR: it appears that recon-all is already running"
    echo "for $subjid based on the presence of $bailfile. It could"
    echo "also be that recon-all was running at one point but"
    echo "died in an unexpected way. If it is the case that there"
    echo "is a process running, you can kill it and start over or"
    echo "just let it run. If the process has died, you should type:"
    echo ""
    echo "rm $bailfile"
    echo ""
    echo "and re-run. Or you can add -no-isrunning to the recon-all"
    echo "command-line. The contents of this file are:"
    echo "----------------------------------------------------------"
    cat  $bailfile
    echo "----------------------------------------------------------"
    exit 1;
  endif
  echo "------------------------------" > $IsRunningFile
  echo "SUBJECT $subjid" >> $IsRunningFile
  echo "HEMI    $hemilist"  >> $IsRunningFile
  echo "DATE `date`"     >> $IsRunningFile
  echo "USER $user"      >> $IsRunningFile
  echo "HOST `hostname`" >> $IsRunningFile
  echo "PROCESSID $$ "   >> $IsRunningFile
  echo "PROCESSOR `uname -m`" >> $IsRunningFile
  echo "OS `uname -s`"       >> $IsRunningFile
  uname -a         >> $IsRunningFile
  echo $VERSION    >> $IsRunningFile
  if($?PBS_JOBID) then
    echo "pbsjob $PBS_JOBID"  >> $IsRunningFile
  endif
endif

# ------- Check FREESURFER_HOME consistency --------------#
set CSDF = $subjdir/scripts/csurfdir
if($DoCleanCSDF) rm -vf $CSDF
if(-e $CSDF) then
  set tmp = `cat $CSDF`;
  if($tmp != $FREESURFER_HOME) then
   echo "INFO: current FREESURFER_HOME does not match that of previous processing." \
     | tee -a $LF
   echo "    Current: $FREESURFER_HOME" | tee -a $LF
   echo "    Previous: $tmp" | tee -a $LF
   sleep 1;
  endif
endif

# --------------- Create the status file ---------------- #
if($#SF == 0) then
  set SF = ($subjdir/scripts/$SF_DEFAULT_NAME)
  if(-e $SF) then
    if(! $AppendStatus) then
       mv -f $SF $SF.old
    else
      echo "\n\n"  >> $SF
      echo "New invocation of recon-all "  >> $SF
      echo "\n\n"  >> $SF
    endif
  endif
else
  if(-e $SF) then
    echo "\n\n"  >> $SF
    echo "New invocation of recon-all "  >> $SF
    echo "\n\n"  >> $SF
  endif
endif
echo "status file for recon-all" >> $SF
date >> $SF

# Put a copy of myself (this script) in the scripts dir
if ("`uname -s`" == "Linux") then
  set FullProgName = `readlink -f $0`; 
  cp -v $FullProgName $subjdir/scripts/recon-all.local-copy
endif


# Wait for a file to appear #
if($#WaitForFile != 0) then
  echo "Waiting for $WaitForFile" |& tee -a $SF |& tee -a $LF
  echo "  WaitSleep $WaitSleep" |& tee -a $SF |& tee -a $LF
  echo "  nWaitsMax $nWaitsMax" |& tee -a $SF |& tee -a $LF
  @ nWaits = 1;
  while(! -e $WaitForFile && $nWaits < $nWaitsMax)
    sleep $WaitSleep;
    @ nWaits = $nWaits + 1;
  end
  if(! -e $WaitForFile ) then
    echo "ERROR: timed out waiting for $WaitForFile"
    goto error_exit;
  endif
  echo "Finished Waiting `date`" |& tee -a $SF |& tee -a $LF
endif

if( ! $RunIt) then
  echo "INFO: -dontrun flag is in effect, so subsequent commands" |& tee -a $LF
  echo "may not have accurate arguments!" |& tee -a $LF
endif


#------------         --------------#
##-----------  -clean --------------#
#------------         --------------#
if($DoCleanSeed) then
  set cmd = ("mv -f $subjdir/scripts/seed-*.crs.man.dat $subjdir/trash");
  echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanCW256) then
  set cmd = (mv -f $subjdir/tmp/cw256 $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanTal) then
  set cmd = (mv -f $subjdir/mri/transforms/talairach.xfm $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
  set cmd = (mv -f $subjdir/mri/orig_nu.mgz $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanLta) then
  set cmd = ("mv -f $subjdir/mri/transforms/*.lta $subjdir/trash")
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanPFH) then
  set cmd = (mv -f $subjdir/mri/optimal_preflood_height $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
  set cmd = (mv -f $subjdir/mri/optimal_skullstrip_invol $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanBM) then
  set cmd = (mv -f $subjdir/mri/brainmask.mgz $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanASeg) then
  set cmd = (mv -f $subjdir/mri/aseg.mgz $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
  set cmd = (mv -f $subjdir/mri/aseg.presurf.mgz $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
  set cmd = (mv -f $subjdir/mri/aseg.manedit.mgz $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanWM) then
  set cmd = (mv -f $subjdir/mri/wm.mgz $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
  set cmd = (mv -f $subjdir/mri/wm.seg.mgz $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanCP) then
  set cmd = (mv -f $subjdir/tmp/control.dat $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanTWM) then
  set cmd = (mv -f $subjdir/tmp/twm.dat $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanBFSE) then
  set cmd = (mv -f $subjdir/mri/brain.finalsurfs.manedit.mgz $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanXopts) then
  set cmd = (mv -f $subjdir/scripts/expert-options $subjdir/trash)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd  |& tee -a $LF
endif
if($DoCleanT2) then
  set mdir = $subjdir/mri
  set flist = ($mdir/T2.mgz $mdir/transforms/T2raw.lta $mdir/transforms/T2raw.auto.lta)
  foreach f ($flist)
    if(-e $f) then
      set cmd = (mv -f $f $subjdir/trash)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $cmd  |& tee -a $LF
    endif
  end
endif
if($DoCleanFLAIR) then
  set mdir = $subjdir/mri
  set flist = ($mdir/FLAIR.mgz $mdir/transforms/FLAIRraw.lta $mdir/transforms/FLAIRraw.auto.lta)
  foreach f ($flist)
    if(-e $f) then
      set cmd = (mv -f $f $subjdir/trash)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $cmd  |& tee -a $LF
    endif
  end
endif


#------------ Handle Seed Points for Fill/Cut, and Watershed ---------------#
set seedfile = $subjdir/scripts/seed-pons.crs.man.dat
if($#PonsSeedCRS) then
  echo "# Manually specified seed CRS for Pons" > $seedfile
  echo $PonsSeedCRS >>  $seedfile
endif
if(-e $seedfile) set PonsSeedCRS = `cat $seedfile | grep -v \#`
set seedfile = $subjdir/scripts/seed-cc.crs.man.dat
if($#CCSeedCRS) then
  echo "# Manually specified seed CRS for CC" > $seedfile
  echo $CCSeedCRS >>  $seedfile
endif
if(-e $seedfile) set CCSeedCRS = `cat $seedfile | grep -v \#`
set seedfile = $subjdir/scripts/seed-lh.crs.man.dat
if($#LHSeedCRS) then
  echo "# Manually specified seed CRS for LH" > $seedfile
  echo $LHSeedCRS >>  $seedfile
endif
if(-e $seedfile) set LHSeedCRS = `cat $seedfile | grep -v \#`
set seedfile = $subjdir/scripts/seed-rh.crs.man.dat
if($#RHSeedCRS) then
  echo "# Manually specified seed CRS for RH" > $seedfile
  echo $RHSeedCRS >>  $seedfile
endif
if(-e $seedfile) set RHSeedCRS = `cat $seedfile | grep -v \#`

set seedfile = $subjdir/scripts/seed-ws.crs.man.dat
if($#WSSeedPoint) then
  echo "# Manually specified seed CRS for watershed" > $seedfile
  echo $WSSeedPoint >>  $seedfile
endif
if(-e $seedfile) set WSSeedPoint = `cat $seedfile | grep -v \#`

#------------ Control Points for Intensity Normalization -----------#
set ControlPointsFile = $subjdir/tmp/control.dat
if(-e $ControlPointsFile) then
  set UseControlPoints = 1;
endif

#------------ Control Points for Registration -----------#
set DefaultTWMControlPointsFile = $subjdir/tmp/twm.dat
if( $UseTWMControlPoints ) then
  #copy twmfile to $subjdir/tmp/twm.dat:
  set cmd = (mkdir -p ${subjdir}/tmp)
  echo "\n $cmd"|& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd |& tee -a $LF
  set cmd = (cp -vf ${TWMControlPointsFile} $DefaultTWMControlPointsFile)
  echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
  if($RunIt) $cmd |& tee -a $LF
  if($status) goto error_exit;
  set TWMControlPointsFile = $DefaultTWMControlPointsFile
else
  if (-e $DefaultTWMControlPointsFile) then
    set UseControlPoints = 1;
    set TWMControlPointsFile = $DefaultTWMControlPointsFile
  endif
endif

#----                      ----#
#---- Conform Width to 256 ----#
##---         -cw256       ----#
if( $DoConformWidth256) then
  echo "-cw256 option is now persistent (remove with -clean-cw256)" \
    |& tee -a $LF
  touch $subjdir/tmp/cw256
endif

#-----------              -----------#
##----------  -show-edits -----------#
#-----------              -----------#
# discover which edits a user has made, list them, and if volume edits were
# made, create a file which shows them (using mri_compile_edits)
if($DoShowEdits) then
  @ edit_count = 0;
  echo "-----------------------------------------------" |& tee -a $LF
  echo "Subject $subjid has the following edits..." |& tee -a $LF

  # control points
  if(-e $ControlPointsFile) then
    @ edit_count = $edit_count + 1;
    set num_cps = (`grep numpoints $ControlPointsFile | awk '{print $2}'`)
    echo "$num_cps control points declared in file $ControlPointsFile" |& tee -a $LF
  endif

  # seeds
  set seedfile = $subjdir/scripts/seed-pons.crs.man.dat
  if(-e $seedfile) then
    @ edit_count = $edit_count + 1;
    echo "Manually specified seed CRS for Pons in file $seedfile" |& tee -a $LF
  endif
  set seedfile = $subjdir/scripts/seed-cc.crs.man.dat
  if($#CCSeedCRS) then
    @ edit_count = $edit_count + 1;
    echo "Manually specified seed CRS for CC in file $seedfile" |& tee -a $LF
  endif
  set seedfile = $subjdir/scripts/seed-lh.crs.man.dat
  if($#LHSeedCRS) then
    @ edit_count = $edit_count + 1;
    echo "Manually specified seed CRS for LH in file $seedfile" |& tee -a $LF
  endif
  set seedfile = $subjdir/scripts/seed-rh.crs.man.dat
  if($#RHSeedCRS) then
    @ edit_count = $edit_count + 1;
    echo "Manually specified seed CRS for RH in file $seedfile" |& tee -a $LF
  endif

  set seedfile = $subjdir/scripts/seed-ws.crs.man.dat
  if($#WSSeedPoint) then
    @ edit_count = $edit_count + 1;
    echo "Manually specified seed CRS for watershed in file $seedfile" \
        |& tee -a $LF
  endif

  # expert opts
  if($#XOptsFile != 0) then
    if (-e $XOptsFile) then
      @ edit_count = $edit_count + 1;
      echo "Expert options declared in file $XOptsFile" |& tee -a $LF
    endif
  endif

  # talairach.xfm
  set xfm = $subjdir/mri/transforms/talairach.xfm
  set xfma = $subjdir/mri/transforms/talairach.auto.xfm
  if(-e $xfm && -e $xfma) then
    diff $xfm $xfma >& /dev/null
    if($status) then
      @ edit_count = $edit_count + 1;
      echo "The talairach.xfm file appears to have been edited" |& tee -a $LF
   endif
  endif

  # cw256
  if(-e $subjdir/tmp/cw256) then
    @ edit_count = $edit_count + 1;
    echo "Conform-width (reduce FOV) to 256 is enabled" |& tee -a $LF
  endif

  # volume edits, as determined by mri_compile_edits:
  # brainmask.mgz, aseg.presurf.mgz, brain.finalsurfs.mgz, wm.mgz, brain.mgz
  set compile_edits=($subjdir/tmp/compile_edits)
  set cmd = (mri_compile_edits ${subjid} $subjdir/mri/edits.mgz)
  if($RunIt) then
    if ( -e $compile_edits) rm -f $compile_edits
    $cmd |& tee -a $LF |& tee -a $compile_edits
    set vol_edits=(`grep mri_compile_edits_found $compile_edits | awk '{print $1}'`)
    if($#vol_edits != 0) then
      @ edit_count = $edit_count + $vol_edits;
    endif
  endif

  # summarize
  echo "$edit_count edits were found for subject $subjid" |& tee -a $LF

endif


#-----------               -----------#
#-----------  Longitudinal -----------#
##---------- -long and -tp ----------#
#-----------               -----------#
if($longitudinal) then
  # if adding a new timepoint was requested by -addtp flag
  if($DoAddTp) then
    echo "Adding $tpNid as timepoint to base $longbaseid" |& tee -a $LF
    set cmd=(mri_add_new_tp $longbaseid $tpNid)
    if($RunIt) $cmd |& tee -a $LF
    if($status) goto error_exit;
  endif

  # init regfile variable with map cross_tp to base:
  # used later in many places
  set tpNtobase_regfile = ${longbasedir}/mri/transforms/${tpNid}_to_${longbaseid}.lta

  # map control.dat file from the cross-sectional data for this subj
  if ( -e ${SUBJECTS_DIR}/${tpNid}/tmp/control.dat && ! $UseLongbaseCtrlVol ) then
    # only if it does not already exist:
    if ( ! -e ${subjdir}/tmp/control.dat ) then
      set cmd = (mkdir -p ${subjdir}/tmp)
      echo "\n $cmd"|& tee -a $LF |& tee -a $CF
      if($RunIt) $cmd |& tee -a $LF
      set cmd = (mri_map_cpdat -in ${SUBJECTS_DIR}/${tpNid}/tmp/control.dat)
      set cmd = ($cmd -out ${subjdir}/tmp/control.dat)
      set cmd = ($cmd -lta $tpNtobase_regfile)
      echo " $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $cmd |& tee -a $LF
      if($status) goto error_exit;
      set UseControlPoints = 1;
    endif
  endif
endif


#-----------                  -----------#
#----------- Convert T1 Input -----------#
##----------    -i <file>     -----------#
#-----------                  -----------#
if($#InputList != 0) then
  @ nth = 1;
  foreach InputVol ($InputList)
    set nthid = `printf %03d.mgz $nth`
    set cmd = (mri_convert $InputVol $subjdir/mri/orig/$nthid)
    $PWD |& tee -a $LF
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
    # sanity-check: make sure each input has the same dimensions
    if ($nth == 1) then
      # assume first input has 'correct' dimensions
      set firstInput = $subjdir/mri/orig/$nthid
      set rows = `mri_info --nrows $firstInput |& tail -n 1`
      set cols = `mri_info --ncols $firstInput |& tail -n 1`
      set slices = `mri_info --nslices $firstInput |& tail -n 1`
    else
      # check subsequent against the first input
      set nextInput = $subjdir/mri/orig/$nthid
      set nrows = `mri_info --nrows $nextInput |& tail -n 1`
      set ncols = `mri_info --ncols $nextInput |& tail -n 1`
      set nslices = `mri_info --nslices $nextInput |& tail -n 1`
      if (($nrows != $rows) || ($ncols != $cols) || ($nslices != $slices)) then
        echo "ERROR: inputs have mismatched dimensions!" |& tee -a $LF
        echo "$firstInput is" |& tee -a $LF
        echo "$rows x $cols x $slices while" |& tee -a $LF
        echo "$InputVol is" |& tee -a $LF
        echo "$nrows x $ncols x $nslices" |& tee -a $LF
        goto error_exit;
      endif
    endif
    @ nth = $nth + 1;
  end # loop over input list
endif

#-----------                          -----------#
#----------- Input T2 or FLAIR images -----------#
##---------- -T2 <file> -FLAIR <file> -----------#
#-----------                          -----------#
if ($DoConvertT2Input || $DoConvertFlairInput) then

  echo "#--------------------------------------------" \
    |& tee -a $LF |& tee -a $CF
  echo "#@# T2/FLAIR Input `date`" \
    |& tee -a $SF |& tee -a $LF |& tee -a $CF

  if ($DoConvertT2Input) then
    set cmd = (mri_convert \
        --no_scale 1 \
        $InputT2Vol \
        $subjdir/mri/orig/T2raw.mgz)
    $PWD |& tee -a $LF
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
  endif

  if ($DoConvertFlairInput) then
    set cmd = (mri_convert \
        --no_scale 1 \
        $InputFlairVol \
        $subjdir/mri/orig/FLAIRraw.mgz)
    $PWD |& tee -a $LF
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
  endif

endif


#-------------                                                 ------------#
#------------- Longitudinal 'base' subject input file creation ------------#
##------------                   -base                         ------------#
#-------------                                                 ------------#
if ($DoCreateBaseInput && $DoCreateBaseSubj) then
  echo "#--------------------------------------------" \
    |& tee -a $LF |& tee -a $CF
  echo "#@# Longitudinal Base Subject Creation `date`" \
    |& tee -a $SF |& tee -a $LF |& tee -a $CF
  if ($#BaseSubjsList == 0) then
    echo "ERROR: must specify subjects for base subject using -base-tp" \
      |& tee -a $LF
    goto error_exit;
  endif
  cd $subjdir > /dev/null
  $PWD |& tee -a $LF
  if($RunIt) rm -f ${BaseSubjsListFname};
  set subjInVols=();
  set normInVols=();
  set ltaXforms=();
  set lta1forms=();
  set ltaAforms=();
  set headInVols=();
  set geodiff=0;
  foreach s ($BaseSubjsList)
    if($RunIt) echo "${s}" >> ${BaseSubjsListFname};
    set normInVols = ($normInVols ${SUBJECTS_DIR}/${s}/mri/norm.mgz)
    set headInVols = ($headInVols ${SUBJECTS_DIR}/${s}/mri/T1.mgz)
    set subjInVols = ($subjInVols ${SUBJECTS_DIR}/${s}/mri/${BaseSubjInvol})
    set ltaname = ${s}_to_${subjid}.lta
    set ltaXforms = ($ltaXforms ${subjdir}/mri/transforms/${ltaname})
    set lta1forms = ($lta1forms ${subjdir}/mri/transforms/${s}_to_${subjid}_norm.lta)
    set ltaAforms = ($ltaAforms ${subjdir}/mri/transforms/${s}_to_${subjid}_affine.lta)

    # check if geometry differs across time
    if ( "$s" != "$BaseSubjsList[1]" ) then
      set cmd = (  mri_diff --notallow-pix --notallow-geo \
                       ${SUBJECTS_DIR}/$s/mri/rawavg.mgz \
                       ${SUBJECTS_DIR}/$BaseSubjsList[1]/mri/rawavg.mgz )
      echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
      if($RunIt) then
        $fs_time $cmd |& tee -a $LF
        if($status) set geodiff = 1;
      endif
    endif
  end

  if ( "$geodiff" == "1" ) then
    echo "\n*******************************************************************************" |& tee -a $LF
    echo "WARNING: Image parameters differ across time, maybe due to aquisition changes?" |& tee -a $LF
    echo "         Consistent changes in, e.g., resolution can potentially bias a " |& tee -a $LF
    echo "         longitudinal study! You can check image parameters by running mri_info" |& tee -a $LF
    echo "         on each input image. Will continue in 10 seconds ..." |& tee -a $LF
    echo "*******************************************************************************\n" |& tee -a $LF
    sleep 10
  endif

  if ($#BaseSubjsList == 1) then
    # if only a single time point, create fake 'base' by making the image upright
    # this assures that also subjects with a single time point get processes as the other
    # subjects in the longitudinal stream

    # 1. make the norm upright (base space)
    set cmd = ( make_upright $normInVols[1] \
        ${subjdir}/mri/norm_template.mgz $ltaXforms[1] )
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif

    # 2. create the upright orig volume
    set cmd = ( mri_convert -rt cubic \
        -at $ltaXforms[1] $subjInVols[1] ${subjdir}/mri/orig.mgz )
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif

  else #more than 1 time point:

    # create the 'mean/median' norm volume:
    set cmd = (mri_robust_template --mov ${normInVols})
    if ($DoAffineBase) then
      # this is only initial rigid reg
      set cmd = ($cmd --lta ${lta1forms})
      set cmd = ($cmd --template ${subjdir}/mri/norm1_template.mgz)
    else
      # this is final rigid reg
      set cmd = ($cmd --lta ${ltaXforms})
      set cmd = ($cmd --template ${subjdir}/mri/norm_template.mgz)
    endif
    set cmd = ($cmd --average ${robust_template_avg_arg})
    set cmd = ($cmd --sat 4.685 )
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif

    if ($DoAffineBase) then
      # use registration on norm.mgz to initialize affine reg on full head imgs:
      set cmd = (mri_robust_template --mov ${headInVols})
      set cmd = ($cmd --lta ${ltaAforms})
      set cmd = ($cmd --average ${robust_template_avg_arg})
      set cmd = ($cmd --template ${subjdir}/mri/head_template.mgz)
      set cmd = ($cmd --sat 4.685 )
      set cmd = ($cmd --ixforms ${lta1forms})
      set cmd = ($cmd --affine)
      echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
      if($RunIt) then
        $fs_time $cmd |& tee -a $LF
        if($status) goto error_exit;
      endif
      # use affine reg to init rigid reg on norm.mgz (fine tuning)
      # not sure if it is really necessary
      set cmd = (mri_robust_template --mov ${normInVols})
      set cmd = ($cmd --lta ${ltaXforms})
      set cmd = ($cmd --average ${robust_template_avg_arg})
      set cmd = ($cmd --template ${subjdir}/mri/norm_template.mgz)
      set cmd = ($cmd --sat 4.685 )
      set cmd = ($cmd --ixforms ${ltaAforms})
      echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
      if($RunIt) then
        $fs_time $cmd |& tee -a $LF
        if($status) goto error_exit;
      endif
    endif

    # create the 'mean/median' input (orig) volume:
    set cmd = (mri_robust_template --mov ${subjInVols})
    set cmd = ($cmd --average ${robust_template_avg_arg})
    set cmd = ($cmd --ixforms ${ltaXforms})
    set cmd = ($cmd --noit)
    set cmd = ($cmd --template ${subjdir}/mri/orig.mgz)
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif

  endif # more than one time point

  # now create the inverse transforms
  cd $subjdir/mri/transforms > /dev/null
  $PWD |& tee -a $LF
  foreach s ($BaseSubjsList)
    set cmd = (mri_concatenate_lta -invert1)
    set cmd = ($cmd ${s}_to_${subjid}.lta)
    set cmd = ($cmd identity.nofile)
    set cmd = ($cmd ${subjid}_to_${s}.lta)
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
  end
  touch $touchdir/base.touch
endif


#-----------             -------------#
##---------- AUTORECON 1 -------------#
#-----------             -------------#


#-----------                            -----------#
#----------- Motion Correct and Average -----------#
##----------        -motioncor          -----------#
#-----------                            -----------#
if($DoMotionCor) then
  echo "#--------------------------------------------" \
    |& tee -a $LF |& tee -a $CF
  echo "#@# MotionCor `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF

  if ($longitudinal) then
    # longitudinal processing to create orig.mgz:

    # in order to create orig.mgz in LONG we need at least 001.mgz in CROSS:
    if ( ! -e ${SUBJECTS_DIR}/${tpNid}/mri/orig/001.mgz ) then
      echo "ERROR: no CROSS run data found in ${SUBJECTS_DIR}/${tpNid}/mri/orig/. Make sure to" \
        |& tee -a $LF
      echo "have a volume called 001.mgz there." |& tee -a $LF
      echo "If you have a second run of data call it 002.mgz, etc." \
        |& tee -a $LF
      echo "See also: http://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/Conversion" \
        |& tee -a $LF
      goto error_exit;
    endif

    # use orig/00?.mgz from cross:
    set CrossList = `ls ${SUBJECTS_DIR}/${tpNid}/mri/orig/[0-9][0-9][0-9].mgz`;
    set origvol = $subjdir/mri/orig.mgz
    set rawvol  = $subjdir/mri/rawavg.mgz

    if($#CrossList == 1) then
      # if only single input, directly resample to base space
      set cmd = (mri_convert -at $tpNtobase_regfile -odt uchar)
      set cmd = ($cmd -rt cubic)
      set cmd = ($cmd ${SUBJECTS_DIR}/${tpNid}/mri/orig/001.mgz)
      set cmd = ($cmd ${origvol})
      echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
      # also create rawavg (usually float) image
      set cmd = (mri_convert -at $tpNtobase_regfile )
      set cmd = ($cmd -rt cubic)
      set cmd = ($cmd ${SUBJECTS_DIR}/${tpNid}/mri/orig/001.mgz)
      set cmd = ($cmd ${rawvol})
      echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
      goto motioncor_post_process
    endif

    # if ltas and iscales exist in cross, copy them over
    set CrossLtas = `ls ${SUBJECTS_DIR}/${tpNid}/mri/orig/[0-9][0-9][0-9].lta`;
    set CrossIscales = `ls ${SUBJECTS_DIR}/${tpNid}/mri/orig/[0-9][0-9][0-9]-iscale.txt`;
    if ( $#CrossLtas > 0 ) then
      # check if one lta for each mgz:
      # here we could better check if really each 00?.mgz has its own 00?.lta in future
      if ( $#CrossList != $#CrossLtas ) then
        echo "ERROR: Orig 00?.mgz runs and number of 00?.lta files must agree in" \
         |& tee -a $LF
        echo "${SUBJECTS_DIR}/${tpNid}/mri/orig/" \
         |& tee -a $LF
        goto error_exit;
      endif
      #copy ltas:
      cd $subjdir/mri/
      set cmd = (cp -vf ${CrossLtas})
      set cmd = ($cmd orig/)
      echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
      #copy iscales if available:
      if ($#CrossIscales > 0) then
        # here we could better check if really each 00?.mgz has its own iscale file in future
        if ( $#CrossList != $#CrossIscales ) then
          echo "ERROR: Orig 00?.mgz runs and number of 00?-iscale.txt files must agree in" \
            |& tee -a $LF
          echo "${SUBJECTS_DIR}/${tpNid}/mri/orig/" \
            |& tee -a $LF
          goto error_exit;
        endif
        cd $subjdir/mri/
        set cmd = (cp -vf ${CrossIscales})
        set cmd = ($cmd orig/)
        echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
        if($RunIt) $fs_time $cmd |& tee -a $LF
        if($status) goto error_exit;
      endif
    else
      # else the ltas don't exist (maybe because 5.0 or fsl was used,
      # they are created by 5.1+ in cross sectional stream)
      # if more than one input, re-receate the ltas:
      if($#CrossList > 1) then
        # set output names for ltas and iscales in long dir
        set LongLtas = ""
        set LongIscales = ""
        foreach nthid ($CrossList)
          set nthname=`basename $nthid .mgz`
          set nthdir=$subjdir/mri/orig
          set nthlta=${nthdir}/${nthname}.lta
          set LongLtas=($LongLtas $nthlta)
          set LongIscales=($LongIscales $nthdir/${nthname}-iscale.txt)
        end
        # perform motion correction again to obtain the ltas (but in long dir, don't touch cross):
        set rawavg  = $subjdir/mri/rawavg.mgz
        # the output rawavg in long will be ignored and not used for anything
        # except maybe debugging, we only need the ltas and iscales!
        set cmd = (mri_robust_template)
        set cmd = ($cmd --mov ${CrossList})
        set cmd = ($cmd --average 1)
        set cmd = ($cmd --template ${rawavg})
        set cmd = ($cmd --satit)
        set cmd = ($cmd --inittp 1)
        set cmd = ($cmd --fixtp)
        set cmd = ($cmd --noit)
        set cmd = ($cmd --iscale)
        set cmd = ($cmd --iscaleout $LongIscales)
        set cmd = ($cmd --subsample 200)
        set cmd = ($cmd --lta $LongLtas)
        echo "#-----------------------------------------------"
        $PWD |& tee -a $LF
        echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
        if($RunIt) then
          $fs_time $cmd |& tee -a $LF
          if($status) goto error_exit;
        endif
        # better get rid of rawavg to avoid confusion
        # as it is not in the base/long space, but in 001.mgz space
        set cmd = (rm -f $rawavg )
        echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
        if($RunIt) then
          $fs_time $cmd |& tee -a $LF
          if($status) goto error_exit;
        endif

      else
        # we should never get here, this case hase been dealt with above
        echo Only single run - should not get here
        goto error_exit;
      endif
    endif

    # now we have the ltas in long (current tp : subjdir)
    set LongLtas = `ls $subjdir/mri/orig/[0-9][0-9][0-9].lta`;
    set LongIscales = `ls $subjdir/mri/orig/[0-9][0-9][0-9]-iscale.txt`;

    # concat ltas (e.g. 002 -> 001 -> base/orig.mgz)
    set ConcatLtas = ""
    foreach nthlta ($LongLtas)
      set nthname=`basename $nthlta .lta`
      set nthdir=$subjdir/mri/orig
      set concatlta=${nthdir}/${nthname}-long.lta
      set ConcatLtas=($ConcatLtas $nthdir/${nthname}-long.lta)
      set cmd = (mri_concatenate_lta $nthlta $tpNtobase_regfile $concatlta)
      echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    end

    # use mri_robust_template just to create the average in base orig space:
    set cmd = (mri_robust_template)
    set cmd = ($cmd --mov ${CrossList})
    set cmd = ($cmd --average 1)
    set cmd = ($cmd --ixforms ${ConcatLtas})
    if ($#LongIscales > 0) then
      set cmd = ($cmd --iscalein ${LongIscales})
    endif
    set cmd = ($cmd --noit)
    set cmd = ($cmd --template ${rawvol})
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif

    # make sure orig is uchar:
    set cmd = (mri_convert -odt uchar ${rawvol} ${origvol})
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif

    goto motioncor_post_process

  endif #longitudinal

  # base processing (we do not need to do anything,
  #   as orig should be there, rawavg is not there, due to difficulties of
  #   averaging different image geometries across time)
  if ($DoCreateBaseSubj) then
    set origvol = $subjdir/mri/orig.mgz
    goto motioncor_post_process
  endif

  # default processing:
  set cmd = ();

  # Get list of input run directories #
  set RunList = ();
  ls $subjdir/mri/orig/[0-9][0-9][0-9].mgz >& /dev/null
  if(! $status) then
    set RunList = `ls $subjdir/mri/orig/[0-9][0-9][0-9].mgz`;
  else
    # No runs found
    ls $subjdir/mri//[0-9][0-9][0-9].mgz >& /dev/null
    if(! $status) then
      set RunList = `ls $subjdir/mri//[0-9][0-9][0-9].mgz`;
    endif
  endif
  if($#RunList == 0 && $RunIt) then
    echo "ERROR: no run data found in $subjdir/mri. Make sure to" \
      |& tee -a $LF
    echo "have a volume called 001.mgz in  $subjdir/mri/orig." |& tee -a $LF
    echo "If you have a second run of data call it 002.mgz, etc." \
      |& tee -a $LF
    echo "See also: http://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/Conversion" \
      |& tee -a $LF
    goto error_exit;
  else
    if ($RunIt) then
      echo "Found $#RunList runs" |& tee -a $LF
      foreach run ($RunList)
        echo $run |& tee -a $LF
      end
    endif
  endif

  # sanity-check: check if the input contains multiple frames, which would
  # be the case if a multi-echo frame mprage was accidently used as input.
  # without this check, andre is unhappy
  foreach RunVol ($RunList)
    echo "Checking for (invalid) multi-frame inputs..." |& tee -a $LF
    set nframes = `mri_info --nframes $RunVol |& tail -n 1`
    if ($nframes != 1) then
      echo "ERROR: input(s) cannot have multiple frames!" |& tee -a $LF
      echo "$RunVol has $nframes frames" |& tee -a $LF
      if ($nframes == 4) then
        echo "If this is a multi-frame MEMPRAGE image," |& tee -a $LF
        echo "use mri_concat --rms to combine echo frames." |& tee -a $LF
      endif
      goto error_exit;
    endif
  end

  set origvol = $subjdir/mri/orig.mgz
  set rawavg  = $subjdir/mri/rawavg.mgz

  # Only one run, copy to rawavg
  if($#RunList == 1) then
    echo "WARNING: only one run found. This is OK, but motion"|& tee -a $LF
    echo "correction cannot be performed on one run, so I'll"|& tee -a $LF
    echo "copy the run to rawavg and continue."|& tee -a $LF
    sleep 2s;
    set cmd = (cp $RunList $rawavg)
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      sleep 1; # Sleep here to assure that they have diff time stamps
      $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
  endif

  # Actually perform the motion correction -- creates rawavg
  if($#RunList > 1) then
    # set names for ltas
    set RunLtas = ""
    set RunIscales = ""
    foreach nthid ($RunList)
      set nthname=`basename $nthid .mgz`
      set nthdir=`dirname $nthid`
      set nthlta=${nthdir}/${nthname}.lta
      set RunLtas=($RunLtas $nthlta)
      set RunIscales=($RunIscales ${nthdir}/${nthname}-iscale.txt)
    end

    if($DoRobustMotionCor) then
      set cmd = (mri_robust_template)
      set cmd = ($cmd --mov ${RunList})
      set cmd = ($cmd --average 1)
      set cmd = ($cmd --template ${rawavg})
      set cmd = ($cmd --satit)
      set cmd = ($cmd --inittp 1)
      set cmd = ($cmd --fixtp)
      set cmd = ($cmd --noit)
      set cmd = ($cmd --iscale)
      set cmd = ($cmd --iscaleout $RunIscales)
      set cmd = ($cmd --subsample 200)
      set cmd = ($cmd --lta $RunLtas)
    else
      # use FSL's flirt to do motion correction (averaging of runs)
      set cmd = (mri_motion_correct.fsl -o $rawavg -wild $RunList)
    endif
    echo "#-----------------------------------------------"
    $PWD |& tee -a $LF
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
      echo $cmd > $touchdir/motion_correct.touch
    endif
  endif

  # At this point, rawavg.mgz exists. Use it to create orig.mgz,
  # conform to COR FOV but keep mgz format.
  set cmd = (mri_convert $rawavg $origvol)
  if($UseCubic) set cmd = ($cmd -rt cubic)
  if($ConformKeepDC) set cmd = ($cmd --conform-dc)
  if($ConformMin) then
    set cmd = ($cmd --conform_min)
  else
    set cmd = ($cmd --conform)
    if($DoConformWidth256 || -e $subjdir/tmp/cw256) then
      # force conform width to 256
      set cmd = ($cmd --cw256)
      # see note just below on when -cw256 is necessary (ie, when FOV > 256)
    endif
  endif
  $PWD |& tee -a $LF
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) then
    $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
    echo $cmd > $touchdir/conform.touch
  endif

  # this is the goto target for longitudinal processing (from above):
  motioncor_post_process:

  # check if FOV > 256 and error exit if so
  set FOV=`mri_info ${origvol} | grep fov: | awk '{print $2}'`
  set FOV_gt_256=`echo "${FOV} > 256" | bc`
  if ($FOV_gt_256) then
    echo "\n****************************************" |& tee -a $LF
    echo "ERROR! FOV=${FOV} > 256" |& tee -a $LF
    echo "Include the flag -cw256 with recon-all!" |& tee -a $LF
    echo "Inspect orig.mgz to ensure the head is fully visible." |& tee -a $LF
    echo "****************************************\n" |& tee -a $LF
    goto error_exit;
  endif

  # Add xfm to orig, even though it does not exist yet. This is a
  # compromise to keep from having to change the time stamp of
  # the orig volume after talairaching.
  set cmd = (mri_add_xform_to_header -c \
             $subjdir/mri/transforms/talairach.xfm \
             $origvol $origvol)
  echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
  if($RunIt) then
    $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
  endif

endif # Motion Correction


#-----------          -----------#
##----------  -deface -----------#
#-----------          -----------#
if($DoDeface) then
  echo "#--------------------------------------------" \
    |& tee -a $LF |& tee -a $CF
  echo "#@# Deface `date`" \
    |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  set cmd = (mri_deface orig.mgz \
     $FREESURFER_HOME/average/$brain_template \
     $FREESURFER_HOME/average/$face_template \
     orig_defaced.mgz)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;
  touch $touchdir/deface.touch
endif


#-----------             -----------#
##----------  -talairach -----------#
#-----------             -----------#
if($DoTalairach) then
talairach:
  echo "#--------------------------------------------" \
    |& tee -a $LF |& tee -a $CF
  echo "#@# Talairach `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  $PWD |& tee -a $LF
  set xfma = transforms/talairach.auto.xfm
  set xfm = transforms/talairach.xfm
  if ($longitudinal) then
    # longitudinal processing:
    if( -e $longbasedir/mri/$xfm) then
      set tal_xfm = $xfm
    else
      set tal_xfm = $xfma
    endif
    # copy from the base (as we are now in same space)
    # if edits were made to base, they will be also in long
    set cmd = (cp $longbasedir/mri/$tal_xfm $subjdir/mri/$xfma)
  else
    # default processing:
    # first run bias correction...
    if ($DoTalairachUseNu) then
      # note that this only works if nu exists
      # this is used as an alternative if tal is giving bad results
      set tal_input = nu.mgz
    else
      # default: run one pass of nu_correct, as its been found that scans
      # with strong bias fields can cause talairach_avi to fail.  note that
      # this orig_nu.mgz uses the full head, whereas nu.mgz is created using
      # the brainmask, and talairach.xfm used to find the white matter ball.
      set tal_input = orig_nu.mgz
      set cmd = (mri_nu_correct.mni --no-rescale --i orig.mgz --o $tal_input)
      # -3T flag support:
      if ($DoNuIntensityCor3T) then
        # 3T params from Zheng, Chee, Zagorodnov 2009 NeuroImage paper
        # "Improvement of brain segmentation accuracy by optimizing
        # non-uniformity correction using N3"
        # namely specifying iterations, proto-iters and distance: 
        set cmd = ($cmd --n 1 --proto-iters 1000 --distance 50)
      else
        # 1.5T default: 
        # DNG: changed distance from 200 to 50 to match 5.3 
        set cmd = ($cmd --n 1 --proto-iters 1000 --distance 50)
      endif
      # allow expert-options override
      set xopts = `fsr-getxopts mri_nu_correct.mni $XOptsFile`;
      set cmd = ($cmd $xopts)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
    # now run the talairach registration...
    if ($UseMincMritotal) then
      # use the MINC mritotal utility
      set xopts = `fsr-getxopts talairach $XOptsFile`;
      set cmd = (talairach --i $tal_input --xfm $xfma $xopts)
    else
      # Avi Snyder's registration tools
      set xopts = `fsr-getxopts talairach_avi $XOptsFile`;
      set cmd = (talairach_avi --i $tal_input --xfm $xfma)
      if($?CustomTalAtlas && "$CustomTalAtlas" != "") then
        # use user-specified atlas found in average dir
        set cmd = ($cmd --atlas "$CustomTalAtlas")
      else
        if($UseYa3tTalAtlas) then
          # special atlas composed of young adults scanned at 3T
          set cmd = ($cmd --atlas 3T18yoSchwartzReactN32_as_orig)
        endif
      endif
      set cmd = ($cmd $xopts)
    endif
  endif
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if ( ! $UseMincMritotal) then
    echo "talairach_avi log file is transforms/talairach_avi.log..." \
        |& tee -a $LF |& tee -a $CF
  endif
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;

  if( -e $xfm && ! $DoCleanTal) then
    echo "\nINFO: $xfm already exists!" \
        |& tee -a $LF |& tee -a $CF
    echo "The new $xfma will not be copied to $xfm" \
        |& tee -a $LF |& tee -a $CF
    echo "This is done to retain any edits made to $xfm" \
        |& tee -a $LF |& tee -a $CF
    echo "Add the -clean-tal flag to recon-all to overwrite $xfm\n" \
        |& tee -a $LF |& tee -a $CF
  endif

  if($DoCleanTal || ! -e $xfm) then
    if(-e $xfm) cp $xfm transforms/bak/talairach.xfm.$DateString
    set cmd = (cp $xfma $xfm)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) then
      sleep 2; # Sleep here to assure that they have diff time stamps
      $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
  endif

  echo $cmd > $touchdir/talairach.touch
endif

# perform the failure detection scheme
if($DoTalCheck) then
  echo "#--------------------------------------------" \
    |& tee -a $LF |& tee -a $CF
  echo "#@# Talairach Failure Detection `date`" \
    |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  # first run Laurence's AFD tool:
  set xopts = `fsr-getxopts talairach_afd $XOptsFile`;
  set xfm = transforms/talairach.xfm
  set cmd = (talairach_afd -T 0.005 -xfm $xfm $xopts)
  set afdTestFailed = 0;
  $PWD |& tee -a $LF
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) then
    $fs_time $cmd |& tee -a $LF
    if($status) then
      set handleTalErrors = 1;
      goto handle_tal_error;
    endif
  endif
  # now run Avi's QA check on the results found in talairach_avi.log
  # see document RLB700_preprocessing_statistics.pdf for details
  set avilog = ($subjdir/mri/transforms/talairach_avi.log)
  if ( ! $UseMincMritotal && -e $avilog) then
    set cmd = (awk -f $FREESURFER_HOME/bin/extract_talairach_avi_QA.awk)
    set cmd = ($cmd $avilog)
    set qalog = ($subjdir/mri/transforms/talairach_avi_QA.log)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    set cmd2 = (tal_QC_AZS $avilog)
    set qalog2 = ($subjdir/mri/transforms/talairach_avi_QA2.log)
    echo "\n $cmd2 \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) then
      rm -f $qalog $qalog2
      $fs_time $cmd >& $qalog
      if($status) then
        echo "ERROR: ${cmd} failed!  See logfile ${qalog}" |& tee -a $LF
        goto error_exit;
      endif
# tal_QC_AZS disabled:
# failing in weird ways. getting caught in infinite 'test' loops.
if (0) then
      $fs_time $cmd2 >& $qalog2
      set tal_QC_AZS_status = $status
      if($tal_QC_AZS_status) then
        echo "INFO: ${cmd2} failed!  See logfile ${qalog2}" |& tee -a $LF
        echo "      skipping tal_QC_AZS test" |& tee -a $LF
        #NJS: skip instead of exit because tal_QC_AZS doesnt work on centos5.4
        #goto error_exit;
      endif
endif
      # first set of QA data:
      set talAviQA=`grep "TalAviQA" ${qalog} | awk '{print $2}'`
      echo "TalAviQA: ${talAviQA}" |& tee -a $LF
      set mean=(0.9781) # buckner40 mean TalAviQA
      set std=(0.0044)  # buckner40 std TalAviQA
      set zscore=`echo "scale=0;(${talAviQA} - ${mean}) / ${std}" | bc`
      echo "z-score: ${zscore}" |& tee -a $LF
      set zscoreThreshold=(-9) # conservative value, but catches M.Harms subjs.
# tal_QC_AZS disabled:
if (0) then
      # second set:
      if ( ! $tal_QC_AZS_status ) then
        set talAviQA2=`grep "atlas_transform_error" ${qalog2} | awk '{print $6}'`
        echo "TalAviQA2: ${talAviQA2}" |& tee -a $LF
        set atlaserrhi=`echo "${talAviQA2} > 24" | bc`
        set atlaserrlo=`echo "${talAviQA2} < -60" | bc`
      else
        set talAviQA2=0
        set atlaserrhi=0
        set atlaserrlo=0
      endif
else
      set talAviQA2=0
      set atlaserrhi=0
      set atlaserrlo=0
endif
      # check scores:
      set handleTalErrors = 0;
      if (($zscore < $zscoreThreshold) || $atlaserrhi || $atlaserrlo) then
        echo "WARNING: Talairach QA check failed!" |& tee -a $LF
        echo "z-score of ${zscore} is <= threshold of ${zscoreThreshold}\n" \
             |& tee -a $LF
# tal_QC_AZS disabled:
if (0) then
        echo "or the atlas xform error of $talAviQA2 is < -60 | > 24" \
             |& tee -a $LF
endif
        set handleTalErrors = 1;
      endif

handle_tal_error:
      if ($handleTalErrors) then
        echo "\nManual Talairach alignment may be necessary, or" |& tee -a $LF
        echo "include the -notal-check flag to skip this test," |& tee -a $LF
        echo "making sure the -notal-check flag follows -all" |& tee -a $LF
        echo "or -autorecon1 in the command string." |& tee -a $LF
        echo "See:\n" |& tee -a $LF
        echo "http://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/Talairach" \
            |& tee -a $LF
        echo "" |& tee -a $LF
        # IF Talairach alignment was run
        # AND it was using the default (711-2C) atlas
        # AND a check failed,
        # THEN first retry avi method using newer 3T Schwartz atlas 
        # THEN (if that fails) try Talairach alignment using MINC mritotal tool
        # ELSE error exit
        if($DoTalairach && ! $UseMincMritotal && ! $longitudinal && \
           ! $UseYa3tTalAtlas ) then
          echo "INFO: Retrying Talairach align using 3T-based atlas...\n" \
            |& tee -a $LF
          set UseYa3tTalAtlas = 1; # try 3T-based atlas
          set UseMincMritotal = 0;
          set DoCleanTal = 1;
          goto talairach;
        else if($DoTalairach && ! $UseMincMritotal && ! $longitudinal) then
          echo "INFO: Trying MINC mritotal to perform Talairach align...\n" \
            |& tee -a $LF
          set UseMincMritotal = 1;
          set DoCleanTal = 1;
          goto talairach;
        else
          echo "\nERROR: Talairach failed!\n" |& tee -a $LF
          goto error_exit;
        endif
      endif
    endif
  endif
endif

#-----------                         -----------#
#----------- Nu Intensity Correction -----------#
##----------     -nuintensitycor     -----------#
#-----------                         -----------#
if($DoNuIntensityCor) then
  echo "#--------------------------------------------" \
    |& tee -a $LF |& tee -a $CF
  echo "#@# Nu Intensity Correction `date`" \
    |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  set xopts = `fsr-getxopts mri_nu_correct.mni $XOptsFile`;
  set cmd = (mri_nu_correct.mni --i orig.mgz --o nu.mgz)
  # to fix a problem with mri_nu_correct.mni messing with the histogramming
  # in some odd subjects (due to voxels outside the brain, which has the
  # effect of altering white matter intensity interpretation downstream),
  # the inclusion of talairach.xfm with mri_nu_correct.mni causes
  # mri_make_uchar to run, which uses the Tal xform to find a ball of voxels
  # that are mostly brain. The top of the intensity histogram in this ball
  # will then be white matter, which allows us to center it at the desired
  # value, approximately (110).
  if ($DoNuMakeUchar) then
    if ( ! -e transforms/talairach.xfm && $RunIt) then
      echo "WARNING: transforms/talairach.xfm does not exist!" \
        |& tee -a $LF
      echo "It will not be used by mri_nu_correct.mni." \
        |& tee -a $LF
    else
      set cmd = ($cmd --uchar transforms/talairach.xfm)
    endif
  endif
  # --cm for hi-res support:
  if ($ConformMin) set cmd = ($cmd --cm)
  # -3T flag support:
  if ($DoNuIntensityCor3T) then
    # 3T params from Zheng, Chee, Zagorodnov 2009 NeuroImage paper
    # "Improvement of brain segmentation accuracy by optimizing
    # non-uniformity correction using N3"
    # namely specifying iterations, proto-iters and distance: 
    set NuIterations = 1;
    set cmd = ($cmd --proto-iters 1000 --distance 50)
  endif
  # add user-specified option (from an 'experts file')
  set cmd = ($cmd --n $NuIterations $xopts)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;
  # Add xfm to nu
  set cmd = (mri_add_xform_to_header -c \
             $subjdir/mri/transforms/talairach.xfm \
             nu.mgz nu.mgz)
  echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
  if($RunIt) then
    $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
  endif
  if($?CREATE_NU_DIFF_FILES) then
    # for debug, create difference volumes between orig and nu files
    set cmd1 = (mri_diff orig.mgz orig_nu.mgz --diff diff-orig-orig_nu.mgz)
    set cmd2 = (mri_diff orig.mgz nu.mgz --diff diff-orig-nu.mgz)
    set cmd3 = (mri_diff orig_nu.mgz nu.mgz --diff diff-orig_nu-nu.mgz)
    echo "\n $cmd1 \n" |& tee -a $LF |& tee -a $CF
    echo "\n $cmd2 \n" |& tee -a $LF |& tee -a $CF
    echo "\n $cmd3 \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) then
      $fs_time $cmd1 |& tee -a $LF
      $fs_time $cmd2 |& tee -a $LF
      $fs_time $cmd3 |& tee -a $LF
    endif
  endif
  # done:
  date > $touchdir/nu.touch
endif

#-----------                          -----------#
#----------- Intensity Normalization1 -----------#
##----------      -normalization      -----------#
#-----------                          -----------#
if($DoNormalization) then
  echo "#--------------------------------------------" \
    |& tee -a $LF |& tee -a $CF
  echo "#@# Intensity Normalization `date`" \
    |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  $PWD |& tee -a $LF
  set xopts = `fsr-getxopts mri_normalize $XOptsFile`;
  if ( $longitudinal && $UseLongbaseCtrlVol) then
    # longitudinal processing stream AND using ctrl_vol.mgz from base:
    if( ! -e $longbasedir/mri/ctrl_vol.mgz || \
        ! -e $longbasedir/mri/bias_vol.mgz ) then
      echo "Recompute intensity norm to create $longbaseid ctrl_vol.mgz" \
        |& tee -a $LF
      set cmd = (mri_normalize)
      if($#Norm3dIters)  set cmd = ($cmd -n $Norm3dIters)
      if($#Norm1_b)      set cmd = ($cmd -b $Norm1_b)
      if($#Norm1_n)      set cmd = ($cmd -n $Norm1_n)
      if($IsMPRAGE)      set cmd = ($cmd -mprage)
      if($IsWashuMPRAGE) set cmd = ($cmd -washu_mprage)
      set cmd = ($cmd \
        -mask $longbasedir/mri/brain.mgz \
        -W $longbasedir/mri/ctrl_vol.mgz $longbasedir/mri/bias_vol.mgz \
        $longbasedir/mri/nu.mgz \
        $longbasedir/mri/T1_tmp.mgz)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
      if( ! -e $longbasedir/mri/ctrl_vol.mgz ) then
        echo "ERROR: unable to generate cntl-point volume for longbase" \
          |& tee -a $LF
        goto error_exit;
      endif
      rm -f $longbasedir/mri/T1_tmp.mgz
    endif
    # use ctrl_vol.mgz from base for normalization:
    set cmd = (mri_normalize \
      -w $subjdir/mri/ctrl_vol.mgz $subjdir/mri/bias_vol.mgz \
      -l $longbasedir/mri/ctrl_vol.mgz $longbasedir/mri/bias_vol.mgz \
      $subjdir/mri/nu.mgz \
      $subjdir/mri/T1.mgz)
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
  else
    if ( $longitudinal ) then
      # the longitudinal stream skips the talairach step, so orig_nu does not
      # get created, so create it now...
      set cmd = (mri_nu_correct.mni --no-rescale --i orig.mgz --o orig_nu.mgz)
      # -3T flag support:
      if ($DoNuIntensityCor3T) then
        # 3T params from Zheng, Chee, Zagorodnov 2009 NeuroImage paper
        # "Improvement of brain segmentation accuracy by optimizing
        # non-uniformity correction using N3"
        # namely specifying iterations, proto-iters and distance: 
        set cmd = ($cmd --n 1 --proto-iters 1000 --distance 50)
      else
      # 1.5T default:
        set cmd = ($cmd --n 1 --proto-iters 1000 --distance 200)
      endif
      # allow expert-options override
      set xopts = `fsr-getxopts mri_nu_correct.mni $XOptsFile`;
      set cmd = ($cmd $xopts)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd >& orig_nu.log
      if($status) goto error_exit;
    endif
    # for cross, base (and long if not useLongBaseCtrlVol) processing streams:
    set cmd = (mri_normalize -g $NormMaxGrad);
    if($UseControlPoints) set cmd = ($cmd -f $ControlPointsFile)
    if($#Norm3dIters)     set cmd = ($cmd -n $Norm3dIters)
    if($#Norm1_b)         set cmd = ($cmd -b $Norm1_b)
    if($#Norm1_n)         set cmd = ($cmd -n $Norm1_n)
    if($IsMPRAGE)         set cmd = ($cmd -mprage)
    if($IsWashuMPRAGE)    set cmd = ($cmd -washu_mprage)
    if($ConformMin)       set cmd = ($cmd -noconform)
    # in base create the ctrl_vol.mgz for init longs later:
    if($DoCreateBaseSubj) set cmd = ($cmd -W ctrl_vol.mgz bias_vol.mgz)
    set cmd = ($cmd $xopts nu.mgz T1.mgz)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
  endif
  echo $cmd > $touchdir/inorm1.touch
endif


#-----------                 -----------#
#----------- Skull Stripping -----------#
##----------   -skullstrip   -----------#
#-----------                 -----------#
if($DoSkullStrip) then
skullstrip:
  echo "#--------------------------------------------" \
    |& tee -a $LF |& tee -a $CF
  echo "#@# Skull Stripping `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  $PWD |& tee -a $LF

  if ($longitudinal) then
    # longitudinal processing stream (copy from base):
    set BM  = brainmask.mgz
    set BMA = brainmask.auto.mgz
    set bmbase = brainmask_${longbaseid}.mgz
    set cmd = (cp -vf ${longbasedir}/mri/${BM} ./${bmbase})
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;

    # first apply mask and keep deletion edits (1)
    set cmd = (mri_mask -keep_mask_deletion_edits)
    set cmd = ($cmd T1.mgz ${bmbase} ${BMA})
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
    # then transfer the 255 edits (still keep deletions)
    set cmd = (mri_mask -transfer 255)
    set cmd = ($cmd -keep_mask_deletion_edits)
    set cmd = ($cmd $BMA ${bmbase} $BMA)
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;

    goto skipped_skullstrip;
  endif

  # base subject stream:
  if($DoCreateBaseSubj) then
    set BMA = brainmask.auto.mgz
    set BM  = brainmask.mgz  #used later
    # create lists with all cross brainmasks and all ltas:
    set bmInVols=""
    set ltaXforms=""
    foreach s ($BaseSubjsList)
      set bmInVols = ($bmInVols ${SUBJECTS_DIR}/${s}/mri/brainmask.mgz)
      set ltaname = ${s}_to_${subjid}.lta
      set ltaXforms = ($ltaXforms ${subjdir}/mri/transforms/${ltaname})
    end

    if ($#bmInVols == 1) then
      # only a single time point in base
      set cmd = (mri_convert \
        -at $ltaXforms[1] \
        -rt nearest \
        $bmInVols[1] brainmask_template.mgz )
      echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
      if($RunIt) then
        $fs_time $cmd |& tee -a $LF
        if($status) goto error_exit;
      endif
    else
      # map all brainmask with nearest neighbor and average 0=mean=logicalOR:
      set cmd = (mri_robust_template --mov ${bmInVols})
      set cmd = ($cmd --average 0)
      set cmd = ($cmd --ixforms ${ltaXforms})
      set cmd = ($cmd --noit)
      set cmd = ($cmd --finalnearest)
      set cmd = ($cmd --template brainmask_template.mgz)
      echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
      if($RunIt) then
        $fs_time $cmd |& tee -a $LF
        if($status) goto error_exit;
      endif
    endif # more than 1 tp

    # create brainmask.auto by applying brainmask_template to T1
    # first apply mask and keep deletion edits (1)
    set cmd = (mri_mask -keep_mask_deletion_edits)
    set cmd = ($cmd T1.mgz brainmask_template.mgz $BMA)
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
    # then copy over the 255 (edits) and still keep deletions (1)
    set cmd = (mri_mask -transfer 255)
    set cmd = ($cmd -keep_mask_deletion_edits)
    set cmd = ($cmd $BMA brainmask_template.mgz $BMA)
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if($RunIt) $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
    goto skipped_skullstrip;
  endif

  # default stream...
  set xopts = `fsr-getxopts mri_watershed $XOptsFile`;
  set opts = ();
  if($WaterShed == 0) set opts = ($opts -n);
  if($WaterShed == 2) set opts = ($opts -wat);
  if($WaterShed == 3) set opts = ($opts -wat+temp);
  if($WaterShed == 4) set opts = ($opts -atlas);
  if($WSLess) set opts = ($opts -less);
  if($WSMore) set opts = ($opts -more);
  if($WSAtlas) set opts = ($opts -atlas);
  if($WSCopy) set opts = ($opts -copy);
  if($#WSSeedPoint != 0) set opts = ($opts -s $WSSeedPoint);
  if($#WSPctPreFlood != 0) then
    set opts = ($opts -h $WSPctPreFlood);
    set DoMultiStrip = 0
  else if( -e optimal_preflood_height) then
    set WSPctPreFlood = `cat optimal_preflood_height`
    set opts = ($opts -h $WSPctPreFlood);
    set DoMultiStrip = 0
    echo "Using optimal preflood height of $WSPctPreFlood" |& tee -a $LF
  endif

  set cmd = (mri_watershed)
  set cmd = ($cmd -rusage $touchdir/rusage.mri_watershed.dat) 
  set BM  = brainmask.mgz

  if ($WSGcaAtlas || $DoMultiStrip) then
    if ( ! $WSUseTalXfm) then # dont bother if using talairach.xfm
      # if using the GCA atlas to help with the skull-strip, then run the
      # GCA registration (mri_em_register) to align to an atlas with a skull,
      # unless that file already exists
      if ( ! -e transforms/talairach_with_skull.lta || $DoCleanLta ) then
        if($NoEMReg == 0) then
          set xopts2 = `fsr-getxopts mri_em_register $XOptsFile`;
          set cmd2 = (mri_em_register)
          if($UseCuda) set cmd2 = (mri_em_register_cuda)
          set cmd2 = ($cmd2 -rusage $touchdir/rusage.mri_em_register.skull.dat) 
          set cmd2 = ($cmd2 -skull)
          set cmd2 = ($cmd2 $xopts2 nu.mgz ${GCADIR}/$GCASkull)
          set cmd2 = ($cmd2 transforms/talairach_with_skull.lta)
        else
          # Convert talairach.xfm to LTA format
          set mni305 = $FREESURFER_HOME/average/mni305.cor.mgz
          set xfm = transforms/talairach.xfm 
          set lta = transforms/talairach_with_skull.lta
          set cmd2 = (lta_convert --src orig.mgz --trg $mni305 \
             --inxfm $xfm  --outlta $lta --subject fsaverage --ltavox2vox)
        endif
        echo "\n $cmd2 \n"|& tee -a $LF |& tee -a $CF
        if($RunIt) $fs_time $cmd2 |& tee -a $LF
        if($status) goto error_exit;
        echo $cmd2 > $touchdir/skull.lta.touch
      endif
    endif

    if ($WSGcaAtlas) then
      # if just using -brain_atlas flag, then include necessary options
      if ($WSUseTalXfm) then
        set opts = (-brain_atlas ${GCADIR}/$GCASkull \
                    transforms/talairach.xfm $opts)
      else
        set opts = (-brain_atlas ${GCADIR}/$GCASkull \
                    transforms/talairach_with_skull.lta $opts)
      endif
    endif
  endif

  if($DoMultiStrip) then

    # when DoMultiStrip is enabled, multiple instances of mri_watershed are
    # run each using the following preflood height parameters
    if ( $?WATERSHED_PREFLOOD_HEIGHTS ) then
      # externally specified in the environment
      set PREFLOOD_HEIGHTS = ( $WATERSHED_PREFLOOD_HEIGHTS )
    else
      # defaults:
      set PREFLOOD_HEIGHTS = ( 5 10 20 30 )
    endif
    set SS_VOLUMES = ( orig orig_nu T1 )
    if ( $?GOTO_LL_CALC ) goto ll_calc

    # create command files, one job for each volume type and preflood height...
    set CMDFS = ()
    foreach vol ( $SS_VOLUMES )
      foreach pfh ( $PREFLOOD_HEIGHTS )
        set cmd = (mri_watershed)
        set cmd = ($cmd -rusage $touchdir/rusage.mri_watershed.dat) 
        set BMA_OLD = brainmask.auto.mgz
        set BMA = brainmask_${vol}_PFH${pfh}.auto.mgz
        if(-e $BMA_OLD && -e $BM && ! $DoCleanBM) then
          set cmd = ($cmd -keep $BMA_OLD $BM $BMA)
        endif
        set cmd = ($cmd $opts -h ${pfh} $xopts ${vol}.mgz $BMA)
        echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
        set CMDF = mri_watershed_${vol}_PFH${pfh}.cmd
        echo "$cmd" > $CMDF
        set CMDFS = ( $CMDFS $CMDF )
      end
    end

    # and launch parallel jobs and wait for completion.
    # the reconbatchjobs script will append each job output to $LF
    if($RunIt) then
      reconbatchjobs $LF $CMDFS
      if($status) then
        # don't exit, since some may encounter 'preflood height' failures
        echo "WARNING: Some multi-skullstrip operations may have failed!" \
            |& tee -a $LF |& tee -a $CF
        echo "         Continuing multi-skullstrip with the others..." \
            |& tee -a $LF |& tee -a $CF
      endif
    endif

    # brainmask volumes actually need to be T1 volumes, because while
    # mri_watershed may run best in some cases using orig.mgz as input,
    # we really want the brainmask to be from T1.mgz (post normalization)
    set CMDFS = ()
    foreach pfh ( $PREFLOOD_HEIGHTS )
      set BMA_ORIG = brainmask_orig_PFH${pfh}.auto.mgz
      if ( -e $BMA_ORIG ) then
        set cmd = (mri_mask T1.mgz $BMA_ORIG $BMA_ORIG)
        echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
        set CMDF = mri_mask_PFH${pfh}.cmd
        echo "$cmd" > $CMDF
        set CMDFS = ( $CMDFS $CMDF )
      endif
    end

    # and launch parallel jobs and wait for completion.
    # the reconbatchjobs script will append each job output to $LF
    if($RunIt) then
      if($#CMDFS != 0) then
        reconbatchjobs $LF $CMDFS
        if($status) goto error_exit;
      endif
    endif
    echo $cmd > $touchdir/skull_strip.touch

    # calculate a measure of skull-strip performance (mri_log_likelihood):
    ll_calc:
    set max_ll = -9999999999;
    set best_pfh = 0;
    set best_vol = ();
    set xopts = `fsr-getxopts mri_log_likelihood $XOptsFile`;
    foreach vol ( $SS_VOLUMES )
      foreach pfh ( $PREFLOOD_HEIGHTS )
        set BMA = brainmask_${vol}_PFH${pfh}.auto.mgz
        if ( ! -e $BMA) continue;
        set LTA = transforms/talairach_with_skull.lta
        set cmd = (mri_log_likelihood -orig T1.mgz $BMA)
        set cmd = ($cmd $xopts ${GCADIR}/$GCA $LTA)
        echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
        if($RunIt) then
          $cmd |& tee -a $LF |& tee -a ll_tmp
          # extract the last line of output, containing the result
          set ll = `tail -n 1 ll_tmp`
          rm -f ll_tmp
        endif
        if($status) goto error_exit;
        # and make the best result the brainmask.auto.mgz
        if($RunIt) then
          echo "$BMA log_likelihood= $ll" |& tee -a $LF
          if ("$ll" == "nan") continue
          if ($ll > $max_ll) then
            set max_ll = $ll;
            set best_pfh = $pfh;
            set best_vol = $vol;
            rm -f brainmask.auto.mgz
            cp $BMA brainmask.auto.mgz
            if($status) goto error_exit;
          endif
        endif
      end
    end
    if($RunIt) then
      if ($best_pfh == 0) then
        echo "ERROR: failure in calculating best preflood height param." \
          |& tee -a $LF
        goto error_exit;
      endif
      echo ""
      echo "Optimal input vol: $best_vol, pre-flood height= $best_pfh, results in log_likelihood= $max_ll" \
        |& tee -a $LF
      #rm -f *_PFH*.auto.*
      #rm -f *_PFH*.*
      # save this result, so that it is used next time, negating the need
      # to run mulitstrip (parallel jobs) again
      echo "$best_vol" > optimal_skullstrip_invol
      echo "$best_pfh" > optimal_preflood_height
    endif
    echo $cmd > $touchdir/log_likelihood.touch
    # make sure this is set properly:
    set BMA = brainmask.auto.mgz

  else

    # single-job skull-strip (default)

    set BMA = brainmask.auto.mgz
    if(-e $BMA && -e $BM && ! $DoCleanBM) set cmd = ($cmd -keep $BMA $BM $BM)
    set INVOL = T1
    if( -e optimal_skullstrip_invol) set INVOL = `cat optimal_skullstrip_invol`
    if("$INVOL" == "T1") set opts = ( -T1 $opts )
    set cmd = ($cmd $opts $xopts $INVOL.mgz $BMA)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
    if("$INVOL" == "orig") then
      # while mri_watershed may run best in some cases using orig.mgz as input,
      # we really want the brainmask to be from T1.mgz (post normalization)
      set cmd = (mri_mask T1.mgz $BMA $BMA)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
    echo $cmd > $touchdir/skull_strip.touch

  endif

  # Gcut:
  # Nanyang Technological University's skull-stripper, which is intended to
  # work best running after watershed.  it further removes dura.  see:
  # 'Skull Stripping Using Graph Cuts', Neuroimage, 2009
  # brainmask.gcuts.mgz is a debug file showing voxels that were cut.
  if($DoGcut) then
    set BMA = brainmask.auto.mgz
    set BGC = brainmask.gcuts.mgz
    set cmd = (rm -f $BGC)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) $cmd |& tee -a $LF
    set xopts = `fsr-getxopts mri_gcut $XOptsFile`;
    set cmd = (mri_gcut $xopts -110 -mult $BMA T1.mgz $BMA $BGC)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    echo "" |& tee -a $LF
    echo "INFO: Care must be taken to thoroughly inspect your data" \
        |& tee -a $LF
    echo "      when using mri_gcut. In particular, inspect the edges of" \
        |& tee -a $LF
    echo "      gm and cerebellum for over-aggressive cutting." |& tee -a $LF
    echo "      Add -segmentation brainmask.gcuts.mgz to the tkmedit" |& tee -a $LF
    echo "      command string to view the voxels which gcut has removed." |& tee -a $LF
    echo "" |& tee -a $LF
    if($RunIt) $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
    # binarize the gcuts file, setting vals to label 999 which appears as red
    # when loaded as a segmentation
    if ( -e $BGC) then
      set cmd = (mri_binarize --i $BGC --o $BGC --binval 999 --min 1)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
  endif

  skipped_skullstrip:

  if( -e $BM && ! $DoCleanBM) then
    echo "\nINFO: brainmask.mgz already exists!" \
        |& tee -a $LF |& tee -a $CF
    echo "The new brainmask.auto.mgz will not be copied to brainmask.mgz." \
        |& tee -a $LF |& tee -a $CF
    echo "This is done to retain any edits made to brainmask.mgz." \
        |& tee -a $LF |& tee -a $CF
    echo "Add the -clean-bm flag to recon-all to overwrite brainmask.mgz.\n" \
        |& tee -a $LF |& tee -a $CF
  endif

  if(! -e $BM || $DoCleanBM) then
    set cmd = (cp $BMA $BM)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) then
      sleep 1; # Sleep here to assure that they have diff time stamps
      $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
  endif

endif



#-----------             -------------#
##---------- AUTORECON 2 -------------#
#-----------             -------------#

#--------------                   ---------------#
#-------------- GCA Registration  ---------------#
##-------------      -gcareg      ---------------#
#--------------                   ---------------#
if($DoGCAReg) then
  echo "#-------------------------------------"|& tee -a $LF |& tee -a $CF
  echo "#@# EM Registration `date`" |& tee -a $SF |& tee -a $LF  |& tee -a $CF
  cd $subjdir/mri > /dev/null
  $PWD |& tee -a $LF
  set xopts = `fsr-getxopts mri_em_register $XOptsFile`;
  if($longitudinal) then
    # longitudinal processing:
     set cmd = (cp -vf $longbasedir/mri/transforms/talairach.lta \
                 transforms/talairach.lta )

  else
    if($NoEMReg == 0) then
      set cmd = (mri_em_register)
      if($UseCuda) set cmd = (mri_em_register_cuda)
      set cmd = ($cmd -rusage $touchdir/rusage.mri_em_register.dat) 
      if($DoCreateBaseSubj) then
        # base subj processing (norm_template instead of nu, and no mask needed):
        set cmd = ($cmd -mask brainmask.mgz $xopts \
                    norm_template.mgz ${GCADIR}/$GCA \
                    transforms/talairach.lta)
      else
        # default processing:
        set lta = transforms/talairach.lta
        if($#GCAOutputName) set lta = transforms/talairach.$GCAOutputName.lta
        set cmd = ($cmd -uns 3 -mask brainmask.mgz $xopts \
                    nu.mgz ${GCADIR}/$GCA $lta)
      endif
    else
      set mni305 = $FREESURFER_HOME/average/mni305.cor.mgz
      set xfm = transforms/talairach.xfm 
      set lta = transforms/talairach.lta
      if($#GCAOutputName) set lta = transforms/talairach.$GCAOutputName.lta
      set cmd = (lta_convert --src orig.mgz --trg $mni305 \
        --inxfm $xfm  --outlta $lta --subject fsaverage --ltavox2vox)
    endif
  endif
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;
  echo $cmd > $touchdir/em_register.touch
endif


#------------                         --------------#
#----------- Canonical Intensity Normalization -----#
##-----------         -canorm         --------------#
#------------                         --------------#
if($DoCANormalize) then
  echo "#--------------------------------------"|& tee -a $LF |& tee -a $CF
  echo "#@# CA Normalize `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  $PWD |& tee -a $LF

  set lta = transforms/talairach.lta
  set norm = norm.mgz
  set ctrl_pts = ctrl_pts.mgz
  if($#GCAOutputName) then
    set lta = transforms/talairach.$GCAOutputName.lta
    set norm = norm.$GCAOutputName.mgz
    set ctrl_pts = ctrl_pts.$GCAOutputName.mgz
  endif

  if ($longitudinal) then
    # longitudinal processing:
    # cp aseg from base to current TP:
    set cmd = (cp -vf ${longbasedir}/mri/aseg.mgz aseg_${longbaseid}.mgz)
    echo "\n $cmd \n" |& tee -a $LF |& tee -a $CF
    if ($RunIt) $fs_time $cmd |& tee -a $LF
    if ($status) goto error_exit;
  endif

  set xopts = `fsr-getxopts mri_ca_normalize $XOptsFile`;
  set cmd = (mri_ca_normalize)
  if($UseCPsWithCaNorm)  set cmd = ($cmd -f $ControlPointsFile)
  if ($longitudinal) then
    # longitudinal processing:
    # use the aseg_base (just created) as initialization of the current TP:
    set cmd = ($cmd -long aseg_${longbaseid}.mgz)
  else
    # default stream
    set cmd = ($cmd -c $ctrl_pts)
  endif
  set cmd = ($cmd -mask brainmask.mgz $xopts)
  if($DoCreateBaseSubj) then
    # base subject stream
    # (use norm vol created during -base-init, no mask needed)
    set cmd = ($cmd norm_template.mgz)
  else
    # default stream
    set cmd = ($cmd nu.mgz)
  endif

  set cmd = ($cmd ${GCADIR}/$GCA $lta $norm)

  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;
  echo $cmd > $touchdir/ca_normalize.touch

endif # DoCANormalize


#------------                        --------------#
#------------ Canonical Registration --------------#
##-----------        -careg          --------------#
#------------                        --------------#
if($DoCAReg) then
  echo "#--------------------------------------"|& tee -a $LF |& tee -a $CF
  echo "#@# CA Reg `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  $PWD |& tee -a $LF

  set xopts = `fsr-getxopts mri_ca_register $XOptsFile`;
  set cmd = (mri_ca_register)
  if($UseCuda) set cmd = (mri_ca_register_cuda)
  set cmd = ($cmd -rusage $touchdir/rusage.mri_ca_register.dat) 
  if($#GCARegIterations) set cmd = ($cmd -n $GCARegIterations)
  if($#GCARegTol) set cmd = ($cmd -tol $GCARegTol)
  set lta = transforms/talairach.lta
  set m3z = transforms/talairach.m3z
  set norm = norm.mgz
  if($#GCAOutputName) then
    set lta = transforms/talairach.$GCAOutputName.lta
    set m3z = transforms/talairach.$GCAOutputName.m3z
    set norm = norm.$GCAOutputName.mgz
  endif
  if($UnCompress) set cmd = ($cmd -uncompress)
  if($BigVentricles) set cmd = ($cmd -bigventricles -smoothness 0.5)
  if( ! $BigVentricles) set cmd = ($cmd -nobigventricles)
  if($DoSecondPassRenorm) set cmd = ($cmd -secondpassrenorm)
  if( ! $longitudinal) set cmd = ($cmd -T $lta)
  if($UseTWMControlPoints) set cmd = ($cmd -twm $TWMControlPointsFile)
  if ( $longitudinal ) then
    # here is tpN's longitudinal command
    # init with warp (m3z) from base
    # no need to concatenate, as we are in same space now
    set cmd = ($cmd -levels 2 -A 1 \
                -l $longbasedir/mri/transforms/talairach.m3z \
                identity.nofile)
  endif

  set cmd = ($cmd $UseCAAlignAfter)
  set cmd = ($cmd -mask brainmask.mgz)
  set cmd = ($cmd $xopts $norm ${GCADIR}/$GCA $m3z)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  set st = $status
  if($st) then
    echo "ERROR: mri_ca_register with non-zero status $st" |& tee -a $LF
    echo "but continuing despite the error" |& tee -a $LF
    #goto error_exit;
  endif

  echo $cmd > $touchdir/ca_register.touch
endif


#------------                                   --------------#
#------------ Inverse of Canonical Registration --------------#
##-----------              -careginv            --------------#
#------------                                   --------------#
if($DoCARegInv) then
  echo "#--------------------------------------"|& tee -a $LF |& tee -a $CF
  echo "#@# CA Reg Inv `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  $PWD |& tee -a $LF
  # set xopts = `fsr-getxopts mri_ca_register $XOptsFile`; # careful here
  set cmd = (mri_ca_register)
  #NJS: cuda version doesnt work will -invert-and-save:
  #if($UseCuda) set cmd = (mri_ca_register_cuda)
  set cmd = ($cmd -invert-and-save transforms/talairach.m3z)
  set cmd = ($cmd -rusage $touchdir/rusage.mri_ca_register.inv.dat) 
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) then
    echo "ERROR: mri_ca_register with non-zero status $status" |& tee -a $LF
    goto error_exit;
  endif
  echo $cmd > $touchdir/ca_register_inv.touch
endif


#------------                                    --------------#
#------------ Removes neck and part of the face  --------------#
##----------- -rmneck, needed if doing SkullLTA  --------------#
#------------                                    --------------#
if($DoRemoveNeck) then
  echo "#--------------------------------------"|& tee -a $LF |& tee -a $CF
  echo "#@# Remove Neck `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  set xopts = `fsr-getxopts mri_remove_neck $XOptsFile`;
  set xform = (transforms/talairach.m3z)
  if ( ! -e $xform) then
    echo "INFO: $xform not found, using talairach.lta instead" \
      |& tee -a $LF |& tee -a $CF
    set xform = (transforms/talairach.lta)
  endif
  set cmd = (mri_remove_neck -radius $RmNeckRadius $xopts \
             nu.mgz $xform \
             ${GCADIR}/$GCA nu_noneck.mgz)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;
  echo $cmd > $touchdir/mri_remove_neck.touch
endif

#------------                                       --------------#
#------------ Recompute lta with skull but no neck  --------------#
##-----------              -skull-lta               --------------#
#------------                                       --------------#
if($DoSkullLTA) then
  echo "#--------------------------------------"|& tee -a $LF |& tee -a $CF
  echo "#@# SkullLTA `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  set xopts = `fsr-getxopts mri_em_register $XOptsFile`;
  set cmd = (mri_em_register)
  if ($UseCuda) set cmd = (mri_em_register_cuda)
  set cmd = ($cmd -skull -t transforms/talairach.lta \
     $xopts nu_noneck.mgz ${GCADIR}/$GCASkull \
     transforms/talairach_with_skull_2.lta)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;
  echo $cmd > $touchdir/skull_2.lta.touch
endif


#--------------                      --------------#
#-------------- SubCort Segmentation --------------#
##-------------        -calabel      --------------#
#--------------                      --------------#
if($DoCALabel) then
  echo "#--------------------------------------"|& tee -a $LF |& tee -a $CF
  echo "#@# SubCort Seg `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  set DoASegMerge = 1; # Force

  if($DoCleanASeg) then
    set cmd = (rm -f aseg.presurf.mgz aseg.manedit.mgz)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $cmd |& tee -a $LF
      if($status) goto error_exit;
  endif

  # ----------- Prior to changing aseg.auto.mgz ------------#
  if( -e aseg.presurf.mgz && -e aseg.auto.mgz ) then
    # aseg.presurf.mgz and aseg.auto.mgz DO exist (at least 2nd pass)
    if( ! -e aseg.manedit.mgz) then
      # aseg.manedit.mgz does NOT exist
      # Check for diffs between aseg.auto and aseg.presurf
      # Note: if there are no diffs, then aseg.manedit.mgz is not created
      set cmd = (mri_seg_diff --seg1 aseg.auto.mgz \
        --seg2 aseg.presurf.mgz --diff aseg.manedit.mgz)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if( ! $RunIt) then
        echo "INFO: mri_seg_diff was not actually executed," |& tee -a $LF
        echo "so subsequent commands (shown with -dontrun)" |& tee -a $LF
        echo "may not have accurate arguments!" |& tee -a $LF
      endif
      if($RunIt) $cmd |& tee -a $LF
      if($status) goto error_exit;
    else
      # aseg.manedit.mgz DOES exist
      # This has to be handled in two stages in case manedit itself has
      # been edited.
      # 1. Check for diffs between aseg.auto and aseg.presurf (new & old edits)
      set cmd = (mri_seg_diff --seg1 aseg.auto.mgz \
        --seg2 aseg.presurf.mgz --diff aseg.manedit.tmp.mgz)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if( ! $RunIt) then
        echo "INFO: mri_seg_diff was not actually executed," |& tee -a $LF
        echo "so subsequent commands (shown with -dontrun)" |& tee -a $LF
        echo "may not have accurate arguments!" |& tee -a $LF
      endif
      if($RunIt) $cmd |& tee -a $LF
      if($status) goto error_exit;
      # 2. Merge new and old edits with manedit. If manedit has not
      # been edited, then there will be no change.
      if(-e aseg.manedit.tmp.mgz) then
        set cmd = (mri_seg_diff --seg aseg.manedit.mgz \
          --diff-in aseg.manedit.tmp.mgz --merged aseg.manedit.mgz)
        echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
        if($RunIt) $cmd |& tee -a $LF
        if($status) goto error_exit;
        rm aseg.manedit.tmp.mgz
      endif
    endif
  endif

  #------------------- mri_ca_label ------------------------#
  # Run mri_ca_label to create aseg.auto_noCCseg.mgz
  set xopts = `fsr-getxopts mri_ca_label $XOptsFile`;
  set rusage = $touchdir/rusage.mri_ca_label.dat 
  set cmd = (mri_ca_label -relabel_unlikely 9 .3 -prior 0.5)
  set cmd = ($cmd $UseCAAlign)
  if($NoWMSA) set cmd = ($cmd -nowmsa)
  set asegbasename=(aseg.auto_noCCseg)
  if($longitudinal && $UseAsegFusion) then
      # longitudinal processing: use 'fused' aseg just created
      set cmd=($cmd -r $subjdir/mri/aseg.fused.mgz)
      set cmd=($cmd -ri $longbasedir/mri/${asegbasename}.label_intensities.txt)
  endif
  set cmd = ($cmd $xopts norm.mgz)
  set cmd = ($cmd transforms/talairach.m3z)
  set cmd = ($cmd ${GCADIR}/$GCA)
  set cmd = ($cmd ${asegbasename}.mgz)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;

  #--- Corpus Callosum Segmentation ---#
  # Run mri_cc to create aseg.auto.mgz
  set xopts = `fsr-getxopts mri_cc $XOptsFile`;
  set cmd = (mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz \
    -lta ${subjdir}/mri/transforms/cc_up.lta $xopts $subjid);
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;

  echo $cmd > $touchdir/ca_label.touch
endif

# --------- Incorporate Manual ASeg Changes ----------------------
# Allow it to be done outside of CALabel. Note that DoASegMerge
# will be 1 if CALabel is run. This just allows the user to
# merge with Man Edits even if CALabel was not run, which will
# only have an effect if the user has actually changed aseg.manedit
if($DoASegMerge) then
  echo "#--------------------------------------"|& tee -a $LF |& tee -a $CF
  echo "#@# Merge ASeg `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  if(-e aseg.manedit.mgz ) then
    set cmd = (mri_seg_diff --seg aseg.auto.mgz \
      --diff-in aseg.manedit.mgz --merged aseg.presurf.mgz)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if( ! $RunIt) then
      echo "INFO: mri_seg_diff was not actually executed," |& tee -a $LF
      echo "so subsequent commands (shown with -dontrun)" |& tee -a $LF
      echo "may not have accurate arguments!" |& tee -a $LF
    endif
    if($RunIt) $cmd |& tee -a $LF
    if($status) goto error_exit;
  else
    set cmd = (cp aseg.auto.mgz aseg.presurf.mgz)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) $cmd |& tee -a $LF
    if($status) goto error_exit;
  endif

  echo $cmd > $touchdir/asegmerge.touch
endif


#-----------                          -----------------#
#----------- Intensity Normalization2 -----------------#
##----------     -normalization2      -----------------#
#-----------                          -----------------#
if($DoNormalization2) then
  echo "#--------------------------------------------" \
    |& tee -a $LF  |& tee -a $CF
  echo "#@# Intensity Normalization2 `date`" \
    |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  set xopts = `fsr-getxopts mri_normalize $XOptsFile`;
  set cmd = (mri_normalize);
  if($UseControlPoints)   set cmd = ($cmd -f $ControlPointsFile)
  if($#Norm3dIters)       set cmd = ($cmd -n $Norm3dIters)
  if($#Norm2_b)           set cmd = ($cmd -b $Norm2_b)
  if($#Norm2_n)           set cmd = ($cmd -n $Norm2_n)
  if($IsMPRAGE)           set cmd = ($cmd -mprage)
  if($IsWashuMPRAGE)      set cmd = ($cmd -washu_mprage)
  if($ConformMin)         set cmd = ($cmd -noconform)
  if($UseAseg && ! $NoAsegInorm2) set cmd = ($cmd -aseg aseg.presurf.mgz)
  if($NoNormMGZ) then
    # norm.mgz doesnt exist with -noaseg or -nosubcortseg, so use brainmask.mgz
    set cmd = ($cmd $xopts brainmask.mgz brain.mgz)
  else
    set cmd = ($cmd -mask brainmask.mgz $xopts norm.mgz brain.mgz)
  endif
  $PWD |& tee -a $LF
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;
  echo $cmd > $touchdir/inorm2.touch
endif


#--------------                        ---------------#
#-------------- Create BrainFinalSurfs ---------------#
##-------------         -maskbfs       ---------------#
#--------------                        ---------------#
if($DoMaskBFS) then
  echo "#--------------------------------------------" \
    |& tee -a $LF  |& tee -a $CF
  echo "#@# Mask BFS `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null
  set xopts = `fsr-getxopts mri_mask $XOptsFile`;
  # Set threshold to WM_MIN_VAL=5
  set cmd = (mri_mask -T 5 $xopts brain.mgz brainmask.mgz brain.finalsurfs.mgz)
  $PWD |& tee -a $LF
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;
  # if they exist, transfer finalsurfs edits (voxels=255 and voxels=1)
  if( -e brain.finalsurfs.manedit.mgz ) then
    set cmd = (mri_mask -transfer 255)
    set cmd = ($cmd -keep_mask_deletion_edits)
    set cmd = ($cmd brain.finalsurfs.mgz brain.finalsurfs.manedit.mgz \
        brain.finalsurfs.mgz)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
  else
    if ( $longitudinal && -e ${SUBJECTS_DIR}/${tpNid}/mri/brain.finalsurfs.manedit.mgz ) then
      # transfer edits from cross
      set cmd = (mri_mask -transfer 255)
      set cmd = ($cmd -keep_mask_deletion_edits)
      set cmd = ($cmd -xform $tpNtobase_regfile)
      set cmd = ($cmd brain.finalsurfs.mgz)
      set cmd = ($cmd ${SUBJECTS_DIR}/${tpNid}/mri/brain.finalsurfs.manedit.mgz)
      set cmd = ($cmd brain.finalsurfs.mgz)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
  endif
endif


#----------------                 -------------------#
#---------------- WM Segmentation -------------------#
##---------------  -segmentation  -------------------#
#----------------                 -------------------#
if($DoSegmentation) then
  echo "#--------------------------------------------" \
    |& tee -a $LF |& tee -a $CF
  echo "#@# WM Segmentation `date`" |& tee -a $SF |& tee -a $LF |& tee -a $CF
  cd $subjdir/mri > /dev/null

  set AllowLongWMtransfers = "1"
  if(! $DoCleanWM && -e wm.mgz) then
    # test if edits were made
    # count 255 edits:
    set cmd = (mri_binarize --i wm.mgz --min 255 --max 255 \
        --o wm255.mgz --count wm255.txt)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) $cmd |& tee -a $LF
    if($status) goto error_exit;
    # count 1 edits:
    set cmd = (mri_binarize --i wm.mgz --min 1 --max 1 \
        --o wm1.mgz --count wm1.txt)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) $cmd |& tee -a $LF
    if($status) goto error_exit;
    # remove unnecessary mgz's
    set cmd = (rm wm1.mgz wm255.mgz)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) $cmd |& tee -a $LF
    if($status) goto error_exit;
    # get number of edits:
    set wm255=""
    if ( -e wm255.txt )  then
      set wm255 = `cat wm255.txt | awk ' { print $1 } '`
    endif
    set wm1=""
    if ( -e wm1.txt ) then
      set wm1 = `cat wm1.txt | awk ' { print $1 } '`
    endif

    # if edits exist:
    if ( "$wm1" != "0" || "$wm255" != "0" ) then
      echo "Found wm edits: $wm1 deletes, $wm255 fills" |& tee -a $LF
      set cmd = (cp wm.mgz wm.seg.mgz)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $cmd |& tee -a $LF
      if($status) goto error_exit;
      # wm.mgz was available (probably edited):
      # don't allow transfers in long below
      set AllowLongWMtransfers = "0"
    endif
  endif
  # ----------- Segment -------------------------
  set xopts = `fsr-getxopts mri_segment $XOptsFile`;
  set cmd = (mri_segment);
  if ($NoThicken) set cmd = ($cmd -thicken 0)
  # note: check for wm.mgz happens so that -dontrun output is correct
  if(! $DoCleanWM && -e wm.mgz && -e wm.seg.mgz) set cmd = ($cmd -keep)
  if($IsMPRAGE)         set cmd = ($cmd -mprage)
  if($IsWashuMPRAGE)    set cmd = ($cmd -washu_mprage)
  set cmd = ($cmd $WMSeg_wlo $WMSeg_ghi $xopts brain.mgz wm.seg.mgz)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;
  # ----------- Edit with ASeg -------------------------
  if($UseAseg) then
    set xopts = `fsr-getxopts mri_edit_wm_with_aseg $XOptsFile`;
    set cmd = (mri_edit_wm_with_aseg)
    if(! $DoCleanWM) set cmd = ($cmd -keep-in)
    set cmd = ($cmd $xopts wm.seg.mgz brain.mgz \
        aseg.presurf.mgz wm.asegedit.mgz)
    echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
    if($RunIt) $fs_time $cmd |& tee -a $LF
    if($status) goto error_exit;
  endif
  # ----------- PreTess -------------------------
  set xopts = `fsr-getxopts mri_pretess $XOptsFile`;
  set wmlabelstring = "wm"
  if($UseAseg) then
    set WM_vol = wm.asegedit.mgz
  else
    set WM_vol = wm.seg.mgz
  endif
  if($NoNormMGZ) then
    set norm_vol = brain.mgz
  else
    set norm_vol = norm.mgz
  endif
  set cmd = (mri_pretess )
  if(! $DoCleanWM && -e wm.mgz) set cmd = ($cmd -keep)
  set cmd = ($cmd $xopts $WM_vol $wmlabelstring $norm_vol wm.mgz)
  echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
  if($RunIt) $fs_time $cmd |& tee -a $LF
  if($status) goto error_exit;

  if($longitudinal && $AllowLongWMtransfers) then
    if ($UseLongbaseWMedits) then
      # transfer wm edits (voxels=255) from longbase wm.mgz to long tpN wm.mgz
      # and -keep_mask_deletion_edits transfers voxel-deletion (voxels=1) edits
      set cmd = (mri_mask -transfer 255)
      set cmd = ($cmd -keep_mask_deletion_edits)
      set cmd = ($cmd wm.mgz $longbasedir/mri/wm.mgz wm.mgz)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    else #default in -long:
      # transfer wm matter edits (voxels=255) from cross wm.mgz to long tpN wm.mgz
      # and -keep_mask_deletion_edits transfers voxel-deletion (voxels=1) edits
      set cmd = (mri_mask -transfer 255)
      set cmd = ($cmd -keep_mask_deletion_edits)
      set cmd = ($cmd -xform $tpNtobase_regfile)
      set cmd = ($cmd wm.mgz ${SUBJECTS_DIR}/${tpNid}/mri/wm.mgz wm.mgz)
      echo "\n $cmd \n"|& tee -a $LF |& tee -a $CF
      if($RunIt) $fs_time $cmd |& tee -a $LF
      if($status) goto error_exit;
    endif
  endif

  echo $cmd > $touchdir/wmsegment.touch
endif
