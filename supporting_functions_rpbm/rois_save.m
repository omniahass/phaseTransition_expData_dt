
clear all
roi = '/Users/gabriellebaxter/Downloads/roi_slice28_rois/roi_slice28_roi_FH1.mat';
% roi = "/Users/gabriellebaxter/Downloads/ro"
load(roi)

dwi_file = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/WCMyofascial1018/FH/derivatives/all/0050ms/dwiec.nii';
dwi = niftiread(dwi_file);
info = niftiinfo(dwi_file);

roi_slice = single(bwroi_rev);
[nx,ny,nz,~] = size(dwi);
roi_array = zeros(nx,ny,nz);
roi_array(:,:,si) = roi_slice;
roi_array = single(roi_array);

info.ImageSize = size(roi_array);
info.PixelDimensions = info.PixelDimensions(1:numel(size(roi_array)));
niftiwrite(roi_array,"FH.nii.gz",info);