warning('This example function will only work on data in the PiP dropbox folder')

%Example on combined dMRI + multi-echo gradient echo data.  
%This works on pip data from the dropbox folder.

%Loads the image and protocol files, then fits MERA, and plots the spectrum.

%The python workbook can be used to nicely plot the saved .mat file.

%full path of the directory containing images, masks, protocol files, etc 
imgpath = '/Users/paddyslator/Dropbox/placentaJhu/t2sdiff/pip0130';
%mask - will fit MERA to the average signal from this
mask = 'placenta_mask';
%load the files from the image directory into matlab structure
t2diff = load_T2MEdiff(imgpath);
%fit peaks with MERA model to the mean signal from the mask
[MERA_out2D, MERA_fitting_out] = MERA_T2MEdiff_fit(t2diff.MERA_data.(mask));
%save the output 
[~, imgdir] = fileparts(imgpath);
%save in the current directory 
%save([imgdir '_' mask '_MERA_fit'],'MERA_out2D');
%contour plot of the output
figure; hold on;
contour(MERA_out2D.T2,...
    MERA_out2D.T,...
    MERA_out2D.S)
ylabel('D (\times 10^{-3} mm^2 s^{-1})')
xlabel('T2^* (s)')        
set(gca,'yscale','log')
        
 
