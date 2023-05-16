function paramindex = GetKernelParameterIndex(modelname,paramname)

strings = GetKernelParameterStrings(modelname);

for i=1:length(strings)
    if strcmp(strings(i), paramname)
        paramindex = i;
        return
    end
end

paramindex = -1;
    


