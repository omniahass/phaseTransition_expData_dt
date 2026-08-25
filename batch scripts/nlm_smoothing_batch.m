clear all

input_folder  = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell/dwi';
output_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell/dwi_smoothed';

% Get all .nii files
files = dir(fullfile(input_folder, '*.nii'));
%%
for i = 226:length(files)
    disp(files(i))
    
    % Full input pat9
    input_file = fullfile(files(i).folder, files(i).name);
    
    % Read NIfTI
    dwi = niftiread(input_file);
    info = niftiinfo(input_file);
    info.Datatype = 'double';
    
    % Apply smoothing
    dwis = dwi_nlm_smoothing(dwi, [], 3, 0.1);


    % dwis = dwi_nlm_smoothing(temp_dwi, [], 3, 0.1);

    disp('FINAL OUTPUT STATS:')
    disp([min(dwis(:)) max(dwis(:))])
    
    % Output path (same filename)
    output_file = fullfile(output_folder, files(i).name);
    
    % Write NIfTI
    niftiwrite(dwis, output_file,info);
    
    fprintf('Processed: %s\n', files(i).name);
end

%%

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

function dwis = dwi_nlm_smoothing_old(dwi,mask,kernel,threshold)
% ----------------------------------------------------------------------- %
% dwi:       input diffusion mri data [nx ny nz ndwi]
% mask:      mask [nx ny nz]
% kernel:    patch size (odd values)
% threshold: if value is [0,1]   represents percentile of voxels in the patch to average (e.g. 0.10->10%) 
%            if value is [1,inf] selects voxels with distance < threshold*min_distance with reference voxel in patch (e.g. 2)  
% ----------------------------------------------------------------------- %
% Ricardo Coronado Leija. November 16th 2021 
% ----------------------------------------------------------------------- %
if(isempty(mask))
mask = ones(size(dwi(:,:,:,1)),'logical');    
end
% from mpnonlocal (Ben's code)
% pad images
if(isscalar(kernel)); kernel = [kernel, kernel, kernel]; end
kernel = kernel + (mod(kernel, 2)-1);  
k      = (kernel-1)/2; 
kx     = k(1); 
ky     = k(2); 
kz     = k(3);
% ksz    = kx*ky*kz;
ksz = prod(kernel);
nend   = round(threshold*ksz); % only used if threshold [0,1]
mask   = padarray(mask       , [kx, ky, kz]   , 'circular');
data   = padarray(single(dwi), [kx, ky, kz, 0], 'circular');
[sx, sy, sz, N] = size(data);
% define a mask that excludes padded values and extract coordinates
mask(1:kx, :, :) = 0;
mask(:, 1:ky, :) = 0;          
mask(:,:, 1:kz) = 0;
mask(sx-kx+1:sx, :, :) = 0;
mask(:, sy-ky+1:sy, :) = 0; 
mask(:, :, sz-kz+1:sz) = 0;
x = []; y = []; z = []; 
for i = kz+1:sz-kz
    [x_, y_] = find(mask(:,:,i) == 1);
    x = [x; x_]; y = [y; y_];  z = [z; i*ones(size(y_))];
end 
x = x(:); y = y(:); z = z(:);
% ======================================================================= %
% prealocate
dwis = zeros(size(data)); 
% start smoothing
Nx = numel(x);
for nn = 1:Nx 
% ======================================================================= %
% from mpnonlocal (Ben's code)    
X  = data(x(nn)-kx:x(nn)+kx,y(nn)-ky:y(nn)+ky,z(nn)-kz:z(nn)+kz,:);
X  = reshape(X, prod(kernel), N);
if nn == 1
    disp('X stats:')
    disp([min(X(:)) max(X(:))])
end
[min_wgs,min_idx] = refine_patch(X, kernel);
% ======================================================================= %
% max distance
wgs_max = min_wgs(end);
% min distance (first is always zero)
wgs_min = min_wgs(2);
if(threshold >= 1)
% remove very bad voxels using the min distance
goodidx = find(min_wgs < threshold*wgs_min); 
elseif(0 < threshold && threshold < 1)
% keep voxels based on percentile
goodidx = 1:nend;
else
error('dwi_nlm_smoothing(): threshold should not be negative')    
end
min_idx = min_idx(goodidx);
min_wgs = min_wgs(goodidx);

if isempty(goodidx)
    error('No good voxels found')
end
if(isscalar(goodidx))
    wgs_sig = X(min_idx,:)';    
else
    % normalizing weights (FIXED)
    wgs_inv = wgs_max - min_wgs;

    if sum(wgs_inv) < 1e-8
        wgs_nrm = ones(size(wgs_inv)) / numel(wgs_inv);
    else
        wgs_nrm = wgs_inv / sum(wgs_inv);
    end

    % compute smoothed signal
    wgs_sig = sum(X(min_idx,:) .* (wgs_nrm * ones(1,N)))';
end% if-else
% ======================================================================= %
if nn == 1
    disp('Sample wgs_sig values:')
    disp(wgs_sig(1:min(5,end)))
end

dwis(x(nn),y(nn),z(nn),:) = wgs_sig;
% ======================================================================= %
end % nn
% ======================================================================= %
dwis  = unpad(dwis,kernel);
end % main
% ======================================================================= %
% from mpnonlocal (Ben's code)    
function [min_wgs,min_idx] = refine_patch_old(data, kernel)
% refval  = data(ceil(prod(kernel)/2),:,:);
center_idx = ceil(size(data,1)/2);
refval = data(center_idx,:);
refval  = repmat(refval,[prod(kernel),1]);
wdists  = 1/size(data,2) * sum((data - refval).^2, 2);
[min_wgs,min_idx] = sort(wdists, 'ascend');
end
function data = unpad(data,kernel)
k = (kernel-1)/2;
data = data(k(1)+1:end-k(1),k(2)+1:end-k(2),k(3)+1:end-k(3),:,:);
end
% ======================================================================= %