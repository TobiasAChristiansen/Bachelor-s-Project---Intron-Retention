#! /usr/bin/bash

#Taking an argument for the input files' directory
new_binders_location=$1
res_location=$2

#Getting all the filenames
SB_I=$(ls "$new_binders_location"/*SB_*NETMHCpan.txt)
WB_I=$(ls "$new_binders_location"/*WB_*NETMHCpan.txt)
SB_II=$(ls "$new_binders_location"/*SB_*NETMHCIIpan.txt)
WB_II=$(ls "$new_binders_location"/*WB_*NETMHCIIpan.txt)

#concatenating
cat $SB_I > $res_location"/"SB_collected_NETMHCpan.txt
cat $WB_I > $res_location"/"WB_collected_NETMHCpan.txt
cat $SB_II > $res_location"/"SB_collected_NETMHCIIpan.txt
cat $WB_II > $res_location"/"WB_collected_NETMHCIIpan.txt
