#!/bin/tcsh -f


set nRuns = 200
#set nRuns = 41

# Directory
set mainDir = /disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/
# path to actual executable code
set exec = hadd

#set outName = newGtag_3BiasTestFixedNonNueCC_T2_ntag00011_expFV_x200
set outName = newGtag_3BiasTestNoNeventsNCCCRatio_ntag00011_expFV_x200

set outDir = /disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/combinedOuts
set combOutName = $outDir/${outName}.root
set inList = ""

# Initialize queue related stuff
set CSH_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/jobFiles
set LOG_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/logFiles
set ERR_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/errFiles
set OUT_DIR=/disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys/withQueue/tempDir/outFiles


foreach nRun (`seq 1 $nRuns`)

set nRunStr = `printf "%04d" $nRun`
set outTempName = $OUT_DIR/${outName}_r${nRunStr}.root
set cshFile = $CSH_DIR/job_r${nRunStr}.sh

#if ( $nRun < 813 || $nRun > 815 ) then
set inList = "$inList $outTempName"
#endif

end

echo $inList
hadd $combOutName $inList
