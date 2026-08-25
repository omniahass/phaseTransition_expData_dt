
% Bayesian RPBM mapping
clear all
close all
clc
addpath('./supporting_functions_rpbm/');
addpath('./supporting_functions_dti/NIFTI/'); % CHANGED FOR WRAPPER: load_untouch_nii lives here

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

% Load dMRI protocol - bvals, dirs
filename = fullfile(parent_folder,'/derivatives/all/0022ms/0022ms.bval');
bvals = load(filename);
filename = fullfile(parent_folder,'/derivatives/all/0022ms/0022ms.bvec');
dirs = load(filename);
% load b0
filename = fullfile(parent_folder,'/derivatives/all/0022ms/dwiec.nii');
nii_b0 = load_untouch_nii(filename); 
b0 = nii_b0.img(:,:,:,1);

% load Dperp
times = {'0022ms','0042ms','0081ms','0156ms','0300ms'};
maps_folder = fullfile(parent_folder,'maps');
nT = length(times);
for t = 1:nT
    if t == 1
        filename = fullfile(maps_folder, ...
            sprintf('%s_RD.nii.gz',times{t}));
    else
        filename = fullfile(maps_folder, ...
            sprintf('%s_RD_withinscanreg.nii.gz',times{t}));
    end
    nii = load_untouch_nii(filename);
    vol = nii.img;
    if t == 1
        [nx, ny, nz] = size(vol);
        Dperp = zeros(nx, ny, nz, nT, 'single');  % preallocate
    end
    Dperp(:,:,:,t) = vol;
end
% load Dpar
times = {'0022ms','0042ms','0081ms','0156ms','0300ms'};
nT = length(times);
for t = 1:nT
    if t == 1
        filename = fullfile(maps_folder, ...
            sprintf('%s_AD.nii.gz',times{t}));
    else
        filename = fullfile(maps_folder, ...
            sprintf('%s_AD_withinscanreg.nii.gz',times{t}));
    end
    nii = load_untouch_nii(filename);
    vol = nii.img;
    if t == 1
        [nx, ny, nz] = size(vol);
        Dpar = zeros(nx, ny, nz, nT, 'single');  % preallocate
    end
    Dpar(:,:,:,t) = vol;
end

% Load sigma
filename = fullfile(parent_folder,'/derivatives/all/magdn/sigma.nii');
nii_sigma = load_untouch_nii(filename); 
info = niftiinfo(filename); 
% Define a mask where you want the RPBM fits
[~, folder_name] = fileparts(parent_folder(1:end-1)); % Remove trailing slash first
filename = fullfile(parent_folder, [folder_name '_reformatdwi.nii.gz']);
nii_mask = load_untouch_nii(filename); 
% mask = true(size(b0));
mask = logical(nii_mask.img);

% Load measurement times [ms] - t
t = [ 22, 42, 81, 156, 300 ];

% Normalize noise level with b0
sigma_b0_norm = nii_sigma.img./(b0+1e-6);

% Estimate noise level in Dpar and Dperp
Dref = [ 1.7 0.9 ];
[ sigma_Dpar , sigma_Dperp ] = RPBMestimation.LLS_noiseprop_DTI_numerical(Dref,bvals,dirs,sigma_b0_norm,mask);

output_folder = fullfile(parent_folder,'rpbm_mapping');
if ~exist(output_folder,"dir") % make map folder if it doesn't exist
    mkdir(output_folder)
end
filename = fullfile(output_folder,'sigma_dperp.nii');
niftiwrite(sigma_Dperp,filename,info);
% % Run RPBM data-driven estimation
% D0_fixed = 2; % Fix D0 for training and estimation, in um^2/ms
% PRoptions.Ntrain = 2e4; % # training data points
% PRoptions.diam_factor = 6.29; % Choose 4 or 6.29 - 4 comes from Novikov Nat Phys
% bounds = [ D0_fixed D0_fixed ; 10 60 ; 0.001 0.05 ]; % Training data distribution: D0 in um^2/ms, a in um, kappa in um/ms
% PRoptions.bounds = bounds; 
% PRoptions.PRdegree = 4; % Degree of polynomial for estimation
% % Input Dperp needs to be 4D [x y z times]
% [ param ] = RPBMestimation.fit4Ddata(Dperp,t,0.8*sigma_Dperp,D0_fixed,mask,PRoptions);
% % [ param ] = RPBMestimation.fit4Ddata(Dperp,t,sigma_Dperp,D0_fixed,mask,PRoptions);

% Run RPBM data-driven estimation
D0_fixed = 2; % Fix D0 for training and estimation, in um^2/ms
PRoptions.Ntrain = 2e4; % # training data points
PRoptions.diam_factor = 6.29; % Choose 4 or 6.29 - 4 comes from Novikov Nat Phys
bounds = [ D0_fixed D0_fixed ; 10 60 ; 0.001 0.05 ]; % Training data distribution: D0 in um^2/ms, a in um, kappa in um/ms
PRoptions.bounds = bounds; 
PRoptions.PRdegree = 4; % Degree of polynomial for estimation

% % Input Dperp needs to be 4D [x y z times] - D0 is a scalar
% [ param ] = RPBMestimation.fit4Ddata(nii_Dperp.img,t,sigma_Dperp,D0_fixed,mask,PRoptions);

% Input Dperp needs to be 4D [x y z times] - D0 is a 3D array (This example is fake but just to test whether this works)
% D0_3D = 2*nii_Dperp.img(:,:,:,1);
% D0_3D = Dpar(:,:,:,end);
D0_3D = mean(Dpar(:,:,:,(end-1):end), 4);
[ param ] = RPBMestimation.fit4Ddata(Dperp,t,1.25*sigma_Dperp,D0_3D,mask,PRoptions);

a_1d = RPBMestimation.vectorize(param.a,mask);
kappa_1d = RPBMestimation.vectorize(param.kappa,mask);
% figure('Position',[2637 686 1676 520])
% subplot(121), histogram(a_1d,linspace(bounds(2,1),bounds(2,2),40)), title 'a', set(gca,'FontSize',25)
% subplot(122), histogram(kappa_1d,linspace(bounds(3,1),bounds(3,2),40)), title '\kappa', set(gca,'FontSize',25)

filename = fullfile(output_folder, 'a.nii');
save_a = param.a;
save_a(~mask) = 0;
niftiwrite(save_a,filename,info);
% niftiwrite(clip(param.a,0,100),filename,info);
% filename = fullfile(output_folder, sprintf('%s_kappa.nii', subject));
% niftiwrite(param.a,filename,info);
filename = fullfile(output_folder, 'kappa.nii');
niftiwrite(param.kappa,filename,info);