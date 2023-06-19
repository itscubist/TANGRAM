#! /bin/tcsh -f

# Change outName and cardName appropriately as desired arguments
set dir = /disk1/disk02/usr6/bbodur/nueoxygenAnalysis/unbinnedLhWithSys
set outDir = /disk2/disk03/usr8/bbodur/nueoxygenAnalysis/unbinnedLhWithSys
#set outName = outputs/testAtmpdSys22a_PrintBinned_withText
#set cardName =  cards/cardAtmpdSysTest.card
#set exec = mainFitSysPrint
#set outName = outputs/newGtagErrorSysTest3_m1sigma
#set cardName =  cards/cardAtmpdSysTest.card
set exec = mainFitTest
set outName = outputs/testMainFit_29sys8t2
set cardName = cards/cardAtmpd29sysTest_20230507.card
#set exec = mainFitReweight
#set outName = outputs/effectAtmFluxratio.root
#set cardName = cards/cardAtmpd1sys_20230419.card
#set exec = mainFitTest

#set exec = mainFitWithScan
#set cardName = cards/cardAtmpd0sys_20230324.card 
#set outName = outputs/sysScan_NoSysTest

# Need root 6 for weighted KDEs
#setenv ROOTSYS /usr/local/sklib_gcc8/root_v6.22.06 
setenv ROOTSYS /disk1/disk02/usr6/bbodur/mySkSoftware/tagged/root6.22.06/root-6.22.06/buildDir 
source $ROOTSYS/bin/thisroot.csh

# Run code
$dir/$exec $dir/$cardName $outDir/$outName 
