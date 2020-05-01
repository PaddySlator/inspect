#!/bin/sh

#  inspect_explore_all.sh
#
#  inspect_explore_all.sh [ options ] img mask
#
#  Created by Paddy Slator on 30/04/2020.

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



inspect_explore_all(){
    #set shortcut to matlab (doesn't matter for cluster version)
    matlab='/Applications/MATLAB_R2019a.app/bin/matlab'
    
    #make sure inspect is on the path and matlab path
    #inspectpath='~/Documents/MATLAB/inspect'
    #cluster version
    inspectpath='/home/pslator/inspect'
    
    export PATH=$inspectpath:$PATH
    export MATLABPATH=$inspectpath
  
    imgs=()
    gradechoinv=()
    masks=()
    kernel=()

    #get the function arguments from the appropriate flags
#    while getopts 'i:g:m:k:' flag; do
#      case "${flag}" in
#        i) imgs=${OPTARG} ;;
#        g) gradechoinv=${OPTARG} ;;
#        m) masks=${OPTARG} ;;
#        k) kernel=${OPTARG}  ;;
#      esac
#    done
    
        
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
    
    for n in {1..3}
    do
        echo $n
    done
   
    for ((n=0; n<=nimg; n++));
      do
        $matlab -nodesktop -r  "try;inspect_explore('${imgs[n]}','$gradechoinv','${masks[n]}','$kernel'); catch; end; quit" > InSpectOutputLogFile
      done

    #cluster version
    #for ((n=0; n<=nimg; n++));
    #do
    #   qsub matlab.multicpufun inspect_explore '${imgs[n]}' '$gradechoinv' '${masks[n]}' '$kernel'
    #done
  
    #echo "$gradechoinv"
    #echo "$masks"
    #echo "$kernel"
    
    
    
# -nojvm -nosplash -nodisplay -nodesktop
#    $matlab -nodesktop -r  "try;inspect_explore('$imgs','$gradechoinv','$masks','$kernel'); catch; end; quit" > InSpectOutputLogFile
    
    
    

    
    
    
    
    

}


inspect_explore_all $*



