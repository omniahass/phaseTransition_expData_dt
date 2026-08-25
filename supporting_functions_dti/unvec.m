function y = unvec(x,mask)
%
% Converts a 2D nm matrix with n the number of measurements per voxel
% and m the number of voxels with mask == true, to a 4D xyzn matrix.
%
% Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
%
dims = [size(mask,1) size(mask,2) size(mask,3)];
if isfloat(x)
    y = NaN([dims(1) dims(2) dims(3) size(x,1)], class(x));
else
    y = zeros([dims(1) dims(2) dims(3) size(x,1)], class(x));
end

for k = 1:size(x,1)
    if isfloat(x)
        Dummy = NaN(dims, class(x));
    else
        Dummy = zeros(dims, class(x));
    end
    Dummy(mask) = x(k,:);
    y(:,:,:,k) = Dummy;
end
end