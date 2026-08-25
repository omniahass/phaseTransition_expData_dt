% DTI mapping using Jelle DTI.m class

clear all
close all
clc

addpath('./supporting_functions_dti/');
diffusion_times = ["0022ms","0042ms","0081ms","0156ms","0300ms"];

% SET PARENT FOLDER HERE, MAKE SURE IT HAS / ON THE END
% CHANGED FOR WRAPPER: run_heal_pipeline.sh sets the HEAL_PARENT environment
% variable, otherwise the hard-coded default below is used
parent_folder = getenv('HEAL_PARENT');
if isempty(parent_folder)
    parent_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/WCMyofascial2075';
end
if ~endsWith(parent_folder,filesep) % this script needs the trailing /
    parent_folder = [parent_folder filesep];
end
fprintf('Parent folder: %s\n',parent_folder);

for i = 1:length(diffusion_times)
    diffusion_time = diffusion_times(i)
    nifti_folder = strcat(parent_folder,'derivatives/all/',diffusion_time,'/');

    if isfolder(nifti_folder)
        
        % Read image
        img = niftiread(strcat(nifti_folder,'dwiec.nii')); 
        img_smoothed = dwi_nlm_smoothing(img, [], 3, 0.1);
        info = niftiinfo(strcat(nifti_folder,'dwiec.nii'));
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        info.Datatype = 'double';
        
        % Read b matrix and bvecs and bvals
        bmat_folder = strcat(parent_folder,'bmatrix/');
        bmatrix = load(strcat(bmat_folder,diffusion_time,'.txt'));
    
        % Remove the b0 at the beginning and end
        img = img_smoothed(:,:,:,2:end-1);
        bmatrix = bmatrix(2:end-1,:);
        
        % mapping using b matrix
        dti = DTI(bmatrix, 'wlls');
        outliers = dti.outliers(img);
        dt = dti.dt(img);
        eig = dti.eig(dt);
        ad = dti.ad(eig);
        rd = dti.rd(eig);
        l2 = dti.l2(eig);
        l3 = dti.l3(eig);
        fa = dti.fa(eig);
    
        maps_folder = strcat(parent_folder,"maps");
        if ~exist(maps_folder,"dir") % make map folder if it doesn't exist
            mkdir(maps_folder)
        end
        niftiwrite(clip(1000*abs(rd),0,3.5),strcat(maps_folder, "/", diffusion_time, "_RD.nii.gz"),info);
        niftiwrite(clip(1000*abs(ad),0,3.5),strcat(maps_folder, "/", diffusion_time, "_AD.nii.gz"),info);
        niftiwrite(clip(1000*abs(l2),0,3.5),strcat(maps_folder, "/", diffusion_time, "_L2.nii.gz"),info);
        niftiwrite(clip(1000*abs(l3),0,3.5),strcat(maps_folder, "/", diffusion_time, "_L3.nii.gz"),info);
        niftiwrite(fa,strcat(maps_folder, "/", diffusion_time, "_FA.nii.gz"),info);
    end
end

function dwis = dwi_nlm_smoothing(dwi, mask, kernel, threshold)

% ----------------------------------------------------------------------- %
% Robust Non-Local Means smoothing (fixed version)
% ----------------------------------------------------------------------- %

if isempty(mask)
    mask = true(size(dwi(:,:,:,1)));    
end

% ensure kernel is 3D odd
if isscalar(kernel)
    kernel = [kernel kernel kernel];
end
kernel = kernel + (mod(kernel,2)-1);

k  = (kernel-1)/2;
kx = k(1); ky = k(2); kz = k(3);

ksz  = prod(kernel);
nend = max(1, round(threshold * ksz));

% pad
mask = padarray(mask, [kx ky kz], 'circular');
data = padarray(double(dwi), [kx ky kz 0], 'circular');

[sx, sy, sz, N] = size(data);

% remove borders
mask(1:kx,:,:) = 0;
mask(:,1:ky,:) = 0;
mask(:,:,1:kz) = 0;
mask(sx-kx+1:sx,:,:) = 0;
mask(:,sy-ky+1:sy,:) = 0;
mask(:,:,sz-kz+1:sz) = 0;

% voxel indices
[x, y, z] = ind2sub(size(mask), find(mask));

% output
dwis = zeros(size(data));

Nx = numel(x);

for nn = 1:Nx

    % extract patch
    X = data(x(nn)-kx:x(nn)+kx, ...
             y(nn)-ky:y(nn)+ky, ...
             z(nn)-kz:z(nn)+kz, :);

    X = reshape(X, [], N);

    % compute distances
    [min_wgs, min_idx] = refine_patch(X);

    % select neighbors
    if threshold >= 1
        wgs_min = min_wgs(2);
        goodidx = find(min_wgs < threshold * wgs_min);

        if isempty(goodidx)
            goodidx = 1;
        end
    else
        goodidx = 1:nend;
    end

    min_idx = min_idx(goodidx);

    % ------------------------------------------------------------------- %
    % ✅ ROBUST SMOOTHING (no fragile weights)
    % ------------------------------------------------------------------- %
    wgs_sig = mean(X(min_idx,:), 1)';

    % assign (safe reshape)
    dwis(x(nn),y(nn),z(nn),:) = reshape(wgs_sig,1,1,1,[]);

end

% unpad
dwis = dwis(1+kx:end-kx, 1+ky:end-ky, 1+kz:end-kz, :);

end

% ----------------------------------------------------------------------- %
% helper
% ----------------------------------------------------------------------- %
function [min_wgs, min_idx] = refine_patch(data)

    center_idx = ceil(size(data,1)/2);
    refval = data(center_idx,:);

    refval = repmat(refval, [size(data,1),1]);

    wdists = mean((data - refval).^2, 2);

    [min_wgs, min_idx] = sort(wdists, 'ascend');

end