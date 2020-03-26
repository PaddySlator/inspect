function [name,ext] = remove_ext_from_nifti(filename)
%Remove the .nii.gz or .nii extension from a nifti file. If the input is
%not a nifti filename then returns the input.

niigz_pos=strfind(filename,'.nii.gz');
nii_pos=strfind(filename,'.nii');

if ~isempty(niigz_pos)
    name=filename(1:(niigz_pos-1));
    ext=filename(niigz_pos:end);
elseif ~isempty(nii_pos)
    name=filename(1:(nii_pos-1));
    ext=filename(nii_pos:end);
else
    name=filename;
    ext=[];
end



end