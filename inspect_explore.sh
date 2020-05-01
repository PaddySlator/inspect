#!/bin/sh

#  inspect_explore.sh
#  
#
#  Created by Paddy Slator on 29/04/2020.

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



inspect_explore(){
    #set shortcut to matlab
    matlab='/Applications/MATLAB_R2019a.app/bin/matlab'
    #make sure inspect is on the path and matlab path
    inspectpath='~/Documents/MATLAB/inspect'
    export PATH=$inspectpath:$PATH
    export MATLABPATH=$inspectpath
  
    imgs=''
    gradechoinv=''
    masks=''
    kernel=''

    #get the function arguments from the appropriate flags
    while getopts 'i:g:m:k:' flag; do
      case "${flag}" in
        i) imgs=${OPTARG} ;;
        g) gradechoinv=${OPTARG} ;;
        m) masks=${OPTARG} ;;
        k) kernel=${OPTARG}  ;;
      esac
    done
    
    
    echo $imgs
    #echo ${OPTARG}
    echo $gradechoinv
    echo $masks
    echo $kernel
        
# -nojvm -nosplash -nodisplay -nodesktop
    $matlab -nodesktop -r  "try;inspect_explore('$imgs','$gradechoinv','$masks','$kernel'); catch; end; quit" > InSpectOutputLogFile
    
    
    

    
    
    
    
    

}


inspect_explore $*

