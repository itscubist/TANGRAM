#! /bin/tcsh -f

# Change outName and cardName appropriately as desired arguments
set dir = /disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys
#set exec = mainFitProf_FreeSys
#set exec = mainFitReweightPMPrint
#set exec = mainFitSysPrint
#set exec = mainFit
#set exec = mainFitWithScanKs2
set exec = mainFitTest

# Need root 6 for weighted KDEs
#setenv ROOTSYS /usr/local/sklib_gcc8/root_v6.22.06 
setenv ROOTSYS /disk1/disk02/usr6/bbodur/mySkSoftware/tagged/root6.22.06/root-6.22.06/buildDir 
source $ROOTSYS/bin/thisroot.csh

# Make code
cd $dir
make clean
make $exec
