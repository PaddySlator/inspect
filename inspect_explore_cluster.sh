#!/bin/sh

#  inspect_explore_all.sh
#
#  inspect_explore_all.sh [ options ] img mask
#
#  Created by Paddy Slator on 01/01/2020.

# Author: Paddy Slator, p.slator@ucl.ac.uk
#
#
# LICENSE
# <inspect toolbox for qMRI analysis>
# Copyright (C) <2020>  <Paddy J. Slator, p.slator@ucl.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.



inspect_explore_cluster(){
    #make sure inspect is on the path and matlab path
    #cluster version
    inspectpath='~/inspect:~/inspect/utilities:~/inspect/examples:~/inspect/subfunctions:~/inspect/plotting:~/inspect/kernels'
    export PATH=$inspectpath:$PATH
    export MATLABPATH=$inspectpath

    #make sure matlab scripts are on the path
    scriptspath='/home/pslator/matlab_scripts/'
    export PATH=$scriptspath:$PATH	
    
  
    imgs=()
    gradechoinv=()
    masks=()
    kernel=()
            
    for n in "$@" ; do
        arg=false
        case "$n" in
            -img) flag='i' ;;
            -grad) flag='g' ;;
            -mask) flag='m' ;;
            -kernel) flag='k' ;;
            *) arg=true ;;
        esac
        #if this is an argument, use the flag that we are on to assign stuff
        if [ $arg = true ] ; then
            case "$flag" in
                i) imgs+=($n) ;;
                g) gradechoinv+=($n) ;;
                m) masks+=($n) ;;
                k) kernel+=($n) ;;
            esac
        fi
    done
        
    for img in ${imgs[@]}
    do
        echo $img
    done
   
    for mask in ${masks[@]}
    do
        echo $mask
    done
       
    nimg="${#imgs[@]}"
    nmask="${#masks[@]}"
    if [nimg ~= nmask]
        then
            echo "Number of images and masks must be the same!" 1>&2
            exit 1
    fi
               	
    
    echo $PATH

    #cluster version
    
    
    #fit to individual scans
#    for ((n=0; n<nimg; n++));
#    do
#        echo ${imgs[n]}
#	echo $gradechoinv
#	echo ${masks[n]}
#	echo $kernel
#	echo inspect_explore
#
#        qsub /home/pslator/matlab_scripts/matlab.multicpufun inspect_explore ${imgs[n]} $gradechoinv ${masks[n]} $kernel
#    done
    
    
    #fit to all scans at once
    #modify file list into the correct input format for inspect_explore.m
    allimgs=$(echo "${imgs[@]}")
    allmasks=$(echo "${masks[@]}")
    #put {" "} around the outside
    allimgs="{\""$allimgs"\"}"
    allmasks="{\""$allmasks"\"}"
    #replace all spaces with ","
    allimgs=$(echo ${allimgs// /\",\"})
    allmasks=$(echo ${allmasks// /\",\"})

    qsub /home/pslator/matlab_scripts/matlab.multicpufun inspect_explore $allimgs $gradechoinv $allmasks $kernel
    
    
  


}


inspect_explore_cluster $*




