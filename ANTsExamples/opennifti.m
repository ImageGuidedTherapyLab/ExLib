
% download nifti tools for matlab
% http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

% load nifti file 
niidata = load_nii('myfile.nii.gz')

% display image
imagesc(niidata.img(:,:,10))

% save image
save_nii(niidata,'newfile.nii.gz')

