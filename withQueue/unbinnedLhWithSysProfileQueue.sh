#!/bin/tcsh -f


set nRuns = 41
# Directory
set mainDir = /disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys
# path to actual executable code
set exec = /disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/mainFitProf_FixedSys


#set outName = atmpdSysProfileQueue1_x41
set outName = atmpdSysProfileQueue_FixedSys_seed813215999_27sysGreaterThan2percent_9subfit_x41
set cardName =  $mainDir/cards/profCards5/prof_


# Initialize queue related stuff
set CSH_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/jobFiles
set LOG_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/logFiles
set ERR_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/errFiles
set OUT_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/outFiles

set QUEUE = atmpd
set JOBLIMIT = 300
set njobs
@ count = 0

foreach nRun (`seq 1 $nRuns`)

set nRunStr = `printf "%04d" $nRun`
set outTempName = $OUT_DIR/${outName}_r${nRunStr}
set cshFile = $CSH_DIR/job_r${nRunStr}_${outName}.sh
set cardTempName = ${cardName}${nRunStr}.card


# create csh file for qsub
echo '#/bin/tcsh -f' >> $cshFile
echo 'setenv ROOTSYS /disk1/disk02/usr6/bbodur/mySkSoftware/tagged/root6.22.06/root-6.22.06/buildDir' >> $cshFile
echo 'source $ROOTSYS/bin/thisroot.csh' >> $cshFile
echo $exec $cardTempName $outTempName >> $cshFile
# end of create csh file for qsub

# Run queue
## limit number of jobs below JOBLIMIT
	 set id = ${nRunStr}_${outName}
	 #set njobs = `qstat  $QUEUE | grep bbodur | wc -l` 
	 set njobs = `pjstat --rg $QUEUE | grep bbodur | wc -l` 
   echo "current jobs: $njobs" 
   while ( $njobs > $JOBLIMIT )
			sleep 1
			#set njobs = `qstat  $QUEUE | grep bbodur | wc -l` 
			set njobs = `pjstat --rg $QUEUE | grep bbodur | wc -l` 
   echo "current jobs: $njobs" 
			echo "current jobs: $njobs" 
   end

# submit a job
    echo qsub -o $LOG_DIR/$id.log -e $ERR_DIR/$id.err -q $QUEUE $cshFile
		#qsub -o $LOG_DIR/$id.log -e $ERR_DIR/$id.err -q $QUEUE $cshFile
		pjsub -o $LOG_DIR/$id.log -e $ERR_DIR/$id.err -L rscgrp=$QUEUE $cshFile
		sleep 1
    @ count ++

endif

end # end of for
