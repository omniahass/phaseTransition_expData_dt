function [y, mask] = vec(x,mask)
%
% Converts a 4D xyzn matrix with n the number of measurements per voxel 
% to a 2D nm matrix with m the number of voxels with mask == true.
%
% Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
%
if ~exist('mask','var')
    mask = ~isnan(x(:,:,:,1));
end
y = zeros([size(x,4) sum(mask(:))], class(x));
for k = 1:size(x,4);
    Dummy = x(:,:,:,k);
    y(k,:) = Dummy(mask(:));
end
end