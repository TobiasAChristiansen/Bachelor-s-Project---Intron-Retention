#! /usr/bin/bash

#Taking an argument to find the correct directory
targetdir=$1

#Removing temporary directories with little to no use after the result is created
rm -rf ${targetdir}01_formatted_data
rm -rf ${targetdir}02_NETMHC_res
rm -rf ${targetdir}03_SB_WB_NETMHC
rm -rf ${targetdir}04_new_binders
rm -rf ${targetdir}05_collected_data