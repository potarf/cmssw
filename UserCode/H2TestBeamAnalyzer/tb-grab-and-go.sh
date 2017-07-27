#!/bin/bash
#rsync run and process for hcaltestbeam2017
#to be run on lxplus7 from CMSSW_8_1_0_pre7/src/UserCode/H2TestBeamAnalyzer directory
#provide run number
#note the workplace and emap


if [ "$#" -ne 1 ]; then
  printf "To run this script, please input testbeam four digit run number.\n"
  printf "  You will need to use the daq password.\n"
  printf "  Example: ./grab-and-go.sh 2921\n"
  exit 1
fi

workplace=/afs/cern.ch/user/c/caleb/public
cmsplace=CMSSW_8_1_0_pre7/src/UserCode/H2TestBeamAnalyzer
emap=EMAP_AC03_26JUL2017_Phase1_RM1RM2_Phase2_RM3_Legacy.txt
run=$1
printf "Rsyncing root file from cmshcaltb02:/data/spool to lxplus7:"$workplace/$cmsplace"\n"
printf "Using emap: $emap\n"
cd $workplace/$cmsplace
rsync -avz -e "ssh -A daq@cmshcaltb02" :/data/spool/HTB_00"$run".root .
./process_run.py -o . -i ./ -r $run -e $emap 

