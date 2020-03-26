function t2diff_struc = load_T2MEdiff(imgdir)




%find the echo time filename
tefile = dir([imgdir '/te_*txt']);
%convert to seconds
te = importdata([imgdir '/' tefile.name]) * 10^-3;
%correct the echo times - T2* decay starts at first echo (which is a spin echo)
%te = te - min(te);


%find the diffusion grads 
gradfile = dir([imgdir '/*protocol*']);
grad = importdata([imgdir '/' gradfile.name]);
grad = grad(:,1:4);

%
gradecho = [repmat(grad,[length(te) 1]) repelem(te,length(grad))];

t2diff_struc.gradecho = gradecho;






%convert gradecho table to protocol
protocol = gradecho2protocol(gradecho);
%
t2diff_struc.protocol = protocol;



%find the magnitude image with all echos filename
imgfile = dir([imgdir '/*abs*alle.nii*']);


%load the nifti with all echo times
t2diffimg = load_untouch_nii([imgdir '/' imgfile.name]);
%make sure it is double
t2diffimg.img = double(t2diffimg.img);

%add to the structure to be returned 
t2diff_struc.t2diffimg = t2diffimg;


%normalise the image to the b=0 and minte volumes
normt2diffimg = t2diffimg;
normt2diffimg.img = normalise_to_b0_and_temin(t2diffimg.img,gradecho);
t2diff_struc.normt2diffimg = normt2diffimg;



%find and load all masks, calculate the average signal in each mask, make a MERA fitting structure for the average signal
maskfile = dir([imgdir '/*mask.nii*']);


for i=1:length(maskfile)
    thismask = load_untouch_nii([imgdir '/' maskfile(i).name]);   
    % make sure it is double
    thismask.img = double(thismask.img);
    
    thismaskname = remove_ext_from_nifti(maskfile(i).name);
    t2diff_struc.masks.(matlab.lang.makeValidName(thismaskname)) = thismask;

   
    
    %average the signal in this mask 
    thismask_av_sig = calculate_mean_signal(t2diffimg.img, thismask.img);
    %transpose to return column vector
    mask_av_sig.(matlab.lang.makeValidName(thismaskname)) = thismask_av_sig';

    % make a data structure for MERA fitting
    for j=1:length(te)
        MERA_data.(matlab.lang.makeValidName(thismaskname)).D(:,j) =  thismask_av_sig(gradecho(:,5) == te(j))';       
    end
    %b-values - rescale
    MERA_data.(matlab.lang.makeValidName(thismaskname)).t = grad(:,4) * 10^-3;
    %echo times
    MERA_data.(matlab.lang.makeValidName(thismaskname)).t2 = te;

    t2diff_struc.MERA_data = MERA_data;

end



t2diff_struc.mask_av_sig = mask_av_sig;



% %load any joint t2-adc parameter maps
% t2diff_file = dir([imgdir '/*VascularBallB0T2_maps*'  ]);
% for i=1:length(t2diff_file)
%     this_t2diff = load_untouch_nii([imgdir '/' t2diff_file(i).name]);
%     this_t2diff_name = remove_ext_from_nifti(t2diff_file(i).name);
%         
%     t2diff_struc.t2adcmaps.(matlab.lang.makeValidName(this_t2diff_name)) = this_t2diff;
% end





end


