%% RPBM fits on example data
clc,clear,close all

root = '/Volumes/labspace-1/Projects/HEAL/rpbm_fitting/2014';
% Load dMRI protocol - bvals, dirs
bvals = load(fullfile(root,'0022ms_dg.bval'));
dirs = load(fullfile(root,'0022ms_dg.bvec'));

nii_b0 = load_untouch_nii(fullfile(root,'b0.nii.gz')); 
nii_Dperp = load_untouch_nii(fullfile(root,'dperp.nii.gz')); 
nii_sigma = load_untouch_nii(fullfile(root,'sigma.nii')); 
% Define a mask where you want the RPBM fits
nii_mask = load_untouch_nii(fullfile(root,'mask.nii')); 
% mask = true(size(nii_b0.img));
mask = logical(nii_mask.img);

% Load measurement times [ms] - t
t = [ 22, 42, 81, 156, 300 ];

% Normalize noise level with b0
sigma_b0_norm = nii_sigma.img./(nii_b0.img+1e-6);

% Estimate noise level in Dpar and Dperp
Dref = [ 1.7 0.9 ];
[ sigma_Dpar , sigma_Dperp ] = RPBMestimation.LLS_noiseprop_DTI_numerical(Dref,bvals,dirs,sigma_b0_norm,mask);

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
D0_3D = 2*nii_Dperp.img(:,:,:,1);
[ param ] = RPBMestimation.fit4Ddata(nii_Dperp.img,t,sigma_Dperp,D0_3D,mask,PRoptions);

a_1d = RPBMestimation.vectorize(param.a,mask);
kappa_1d = RPBMestimation.vectorize(param.kappa,mask);
figure('Position',[2637 686 1676 520])
subplot(121), histogram(a_1d,linspace(bounds(2,1),bounds(2,2),40)), title 'a', set(gca,'FontSize',25)
subplot(122), histogram(kappa_1d,linspace(bounds(3,1),bounds(3,2),40)), title '\kappa', set(gca,'FontSize',25)

IMGUI(cat(4,mask*100,param.a),[10 80])
IMGUI(cat(4,mask*100,param.kappa),[0 0.05])





%% Example code to run RPBM fit on synthetic data
%% Running RPBM noise propagationon simulated 4D data
clc,clear,close all

D0_gt = 2.6;
% Define measurements [ms]
% t = [ 50 100 200 400 ];
t = [ 22, 42, 81, 156, 300 ];
% Compute training data
Ntest = 5e3; diam_factor = 4;
M = Ntest;
bounds = [ D0_gt D0_gt ; 10 50 ; 0.001 0.05 ]; % D0 in um^2/ms, a in um, kappa in um/ms
[D0_tr,a_tr,kappa_tr,tau_tr,zeta_tr] = RPBMestimation.sample_rpbm_param(bounds,M,diam_factor);
tnorm_all = t./tau_tr;
zeta_all = repmat(zeta_tr,1,length(t));
Dperp_all = D0_tr.*RPBMestimation.get_Dt_RPBM_vectorized(tnorm_all, zeta_all);

SNR = 20; Dperp_ref = 0.8; D0_fixed = D0_gt; mask = [];
sigma_Dperp = Dperp_ref/SNR;
Dperp_meas = Dperp_all + sigma_Dperp*randn(Ntest,length(t));

PRoptions.Ntrain = 1e4;
PRoptions.diam_factor = 4;
PRoptions.bounds = bounds;
PRoptions.PRdegree = 4;
% Input needs to be 4D
Dperp_meas = reshape(Dperp_meas,[10 10 Ntest/100 length(t)]);
tic
[ param ] = RPBMestimation.fit4Ddata(Dperp_meas,t,sigma_Dperp,D0_fixed,mask,PRoptions);
toc
figure('Position',[2686 690 1608 458])
subplot(131), plot(a_tr,param.a(:),'.',bounds(2,:),bounds(2,:),'r-.'), set(gca,'FontSize',25)
xlabel('gt a'), ylabel('pr fits a'), set(gca,'FontSize',25), title(['SNR = ',num2str(SNR),'']), axis([bounds(2,:), bounds(2,:)])
subplot(132), plot(kappa_tr,param.kappa(:),'.',bounds(3,:),bounds(3,:),'r-.'), set(gca,'FontSize',25)
xlabel('gt \kappa'), ylabel('pr fits \kappa'), set(gca,'FontSize',25), title(['SNR = ',num2str(SNR),'']), axis([bounds(3,:), bounds(3,:)])
subplot(133), plot(param.a(:),param.kappa(:),'.',bounds(2:3,1),bounds(2:3,2),'r-.'), set(gca,'FontSize',25)
xlabel('pr fits a'), ylabel('pr fits \kappa'), set(gca,'FontSize',25), title(['SNR = ',num2str(SNR),'']), axis([bounds(2,:), bounds(3,:)])


%% Looking at the noise map
clc,clear,close all

root = '/Volumes/labspace-1/Projects/HEAL/rpbm_fitting/2014';
% Load dMRI protocol - bvals, dirs
bvals = load(fullfile(root,'0022ms_dg.bval'));
dirs = load(fullfile(root,'0022ms_dg.bvec'));

nii_b0 = load_untouch_nii(fullfile(root,'b0.nii.gz')); 
nii_Dperp = load_untouch_nii(fullfile(root,'dperp.nii.gz')); 
nii_sigma = load_untouch_nii(fullfile(root,'sigma.nii')); 
% Define a mask where you want the RPBM fits
nii_mask = load_untouch_nii(fullfile(root,'mask.nii')); 
% mask = true(size(nii_b0.img));
mask = logical(nii_mask.img);

nii_dwi = load_untouch_nii(fullfile(root,'0022ms_dwi.nii'));
b0 = nii_dwi.img(:,:,:,bvals<50);
sigma_diff_3D = sqrt(pi)/2*abs(b0(:,:,:,1)-b0(:,:,:,2));
sigma_diff_3D_sm = imgaussfilt3(sigma_diff_3D, 1.5);
% sigma_diff = imgaussfilt3(sigma_diff, 1.5);
% IMGUI(cat(4,mask,nii_sigma.img./sigma_diff_3D),[0 2])
% IMGUI(cat(4,mask,nii_sigma.img,sigma_diff_3D_sm,sigma_diff_3D),[0 20])

sigma_MP = RPBMestimation.vectorize(nii_sigma.img./mean(b0,4),mask);
sigma_diff = RPBMestimation.vectorize(sigma_diff_3D./mean(b0,4),mask);
sigma_diff_sm = RPBMestimation.vectorize(imgaussfilt3(sigma_diff_3D, 1.5)./mean(b0,4),mask);

bins = linspace(0,0.4,80);
figure('Position',[2637 763 1663 443]), hold on, set(gca,'FontSize',25)
histogram(sigma_MP,bins,'DisplayStyle','Stairs','LineWidth',2)
histogram(sigma_diff,bins,'DisplayStyle','Stairs','LineWidth',2)
histogram(sigma_diff_sm,bins,'DisplayStyle','Stairs','LineWidth',2)
legend('MPPCA','b_0 repeats','b_0 repeats (smoothed)'), title '\sigma noise DWI estimates'


