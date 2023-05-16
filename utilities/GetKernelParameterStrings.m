function strings = GetKernelParameterStrings(kernelname)


switch kernelname
    case 'DT2'
        strings = {'d','t2'};
    case 'DT2T1'
        strings = {'d','t2','t1'};
    case 'ZeppelinT2T1'
        strings = {'dpar','dperp','theta','phi','t2','t1'};
    case 'DD'
        strings = {'d1','d2'};
    case 'D'
        strings = {'d'};
    case 'DK'
        strings = {'d','k'};
    case 'DKT2T1'
        strings = {'d','k','t2','t1'};
    case 'T2'
        strings = {'t2'};
    case 'PAD'
        strings = {'dpar','dperp'};
    case 'PADT2'
        strings = {'dpar','dperp','t2'};
    case 'T2T1'
        strings = {'t2','t1'};
end
        
                      

