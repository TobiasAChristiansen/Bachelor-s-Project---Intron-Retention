#! /usr/bin/bash

#Taking an argument for the input files' directories
NETMHC_location=$1
res_location=$2


#Identifying the NETMHCpan result files and grabbing the SB/WB data
NETMHCpan_files=`ls $NETMHC_location`
for file in ${NETMHCpan_files[@]}; do
    grep "<= SB" $NETMHC_location/$file > $res_location/SB_$file
    grep "<= WB" $NETMHC_location/$file > $res_location/WB_$file
done


#Identifying the NETMHCIIpan result files and grabbing the SB/WB data
#NETMHCIIpan_files=`ls $NETMHC_location | grep *NETMHCIIpan*`
#for file in ${NETMHCIIpan_files[@]}; do
#    grep "<= SB" $NETMHC_location/$file > $res_location/SB_$file
#    grep "<= WB" $NETMHC_location/$file > $res_location/WB_$file
#done