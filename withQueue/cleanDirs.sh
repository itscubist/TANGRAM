#!/bin/tcsh -f

# Delete queue related stuff when the job is finished
set CSH_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/jobFiles
set LOG_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/logFiles
set ERR_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/errFiles
set OUT_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/outFiles

/bin/rm $CSH_DIR/*
/bin/rm $LOG_DIR/*
/bin/rm $ERR_DIR/*
#/bin/rm $OUT_DIR/*
