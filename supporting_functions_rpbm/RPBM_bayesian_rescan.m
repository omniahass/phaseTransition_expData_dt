%% RPBM fits on example data
clc,clear,close all

parent_directory = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell_rescan/';
nifti_folder = fullfile(parent_directory,'dwi_smoothed');
bvecbval_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell/bvecbval';
maps_folder = fullfile(parent_directory,'maps_smoothed');
sigma_folder = fullfile(parent_directory,'sigma');
mask_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell/roi_resampled';
output_folder = fullfile(parent_directory,'rpbm_maps_d0_156and300ms');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

subjects = dir(nifti_folder);
subjects = {subjects.name};
% Cornell
subjects = subjects(contains(subjects, 'WCMyofascial20'));
% subjects = subjects(~contains(subjects, 'WCMyofascial2067'));
baseNames = cellfun(@(x) split(x, '_'), subjects, 'UniformOutput', false);
baseNames = cellfun(@(x) x{1}, baseNames, 'UniformOutput', false);
uniqueNames = unique(baseNames);
temp = uniqueNames';

% % NYU
% subjects = subjects(contains(subjects, 'C2'));
% baseNames = cellfun(@(x) split(x, '_'), subjects, 'UniformOutput', false);
% baseNames = cellfun(@(x) strjoin(x(1:2), '_'), baseNames, 'UniformOutput', false);
% uniqueNames = unique(baseNames);
% mask_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/nyu/roi_resampled';
% % temp = uniqueNames';


%%
for i = 31:length(uniqueNames)
    subject = uniqueNames{i};
    disp(subject)
    
    % Check if ROI exists
    filename = fullfile(mask_folder, sprintf('%s.nii.gz', subject));
    if exist(filename,'file')

    
        % Load dMRI protocol - bvals, dirs
        filename = fullfile(bvecbval_folder, sprintf('%s_0022ms.bval', subject));
        bvals = load(filename);
        filename = fullfile(bvecbval_folder, sprintf('%s_0022ms.bvec', subject));
        dirs = load(filename);
        
        % load b0
        filename = fullfile(nifti_folder, sprintf('%s_0022ms.nii', subject));
        nii_b0 = load_untouch_nii(filename); 
        b0 = nii_b0.img(:,:,:,1);
        % load Dperp
        times = {'0022ms','0042ms','0081ms','0156ms','0300ms'};
        nT = length(times);
        for t = 1:nT
            filename = fullfile(maps_folder,sprintf('%s_%s_RD_reg2wcmscan1.nii.gz', subject, times{t}));
            disp(filename)
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
            filename = fullfile(maps_folder,sprintf('%s_%s_AD_reg2wcmscan1.nii.gz', subject, times{t}));
            nii = load_untouch_nii(filename);
            vol = nii.img;
            if t == 1
                [nx, ny, nz] = size(vol);
                Dpar = zeros(nx, ny, nz, nT, 'single');  % preallocate
            end
            Dpar(:,:,:,t) = vol;
        end
        % Load sigma
        filename = fullfile(sigma_folder, sprintf('%s_sigma_reg2wcmscan1.nii', subject));
        nii_sigma = load_untouch_nii(filename); 
        info = niftiinfo(filename); 
        % Define a mask where you want the RPBM fits
        filename = fullfile(mask_folder, sprintf('%s.nii.gz', subject));
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
        % 
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
        % 
        % % Input Dperp needs to be 4D [x y z times] - D0 is a scalar
        % [ param ] = RPBMestimation.fit4Ddata(nii_Dperp.img,t,sigma_Dperp,D0_fixed,mask,PRoptions);

        % Input Dperp needs to be 4D [x y z times] - D0 is a 3D array (This example is fake but just to test whether this works)
        % D0_3D = 2*nii_Dperp.img(:,:,:,1);
        % D0_3D = Dpar(:,:,:,end);
        D0_3D = mean(Dpar(:,:,:,(end-1):end), 4);
        [ param ] = RPBMestimation.fit4Ddata(Dperp,t,sigma_Dperp,D0_3D,mask,PRoptions);

        % a_1d = RPBMestimation.vectorize(param.a,mask);
        % kappa_1d = RPBMestimation.vectorize(param.kappa,mask);
        % % figure('Position',[2637 686 1676 520])
        % % subplot(121), histogram(a_1d,linspace(bounds(2,1),bounds(2,2),40)), title 'a', set(gca,'FontSize',25)
        % subplot(122), histogram(kappa_1d,linspace(bounds(3,1),bounds(3,2),40)), title '\kappa', set(gca,'FontSize',25)

    
        filename = fullfile(output_folder, sprintf('%s_a.nii', subject));
        save_a = param.a;
        save_a(~mask) = 0;
        niftiwrite(save_a,filename,info);
        % niftiwrite(clip(param.a,0,100),filename,info);
        % filename = fullfile(output_folder, sprintf('%s_kappa.nii', subject));
        % niftiwrite(param.a,filename,info);
        filename = fullfile(output_folder, sprintf('%s_kappa.nii', subject));
        niftiwrite(param.kappa,filename,info);
        
        % 
        % a_1d = RPBMestimation.vectorize(param.a,mask);
        % kappa_1d = RPBMestimation.vectorize(param.kappa,mask);
        % figure('Position',[2637 686 1676 520])
        % subplot(121), histogram(a_1d,linspace(bounds(2,1),bounds(2,2),40)), title 'a', set(gca,'FontSize',25)
        % subplot(122), histogram(kappa_1d,linspace(bounds(3,1),bounds(3,2),40)), title '\kappa', set(gca,'FontSize',25)
        % 
        % IMGUI(cat(4,mask*100,param.a),[10 80])
        % IMGUI(cat(4,mask*100,param.kappa),[0 0.05])

    end
end


%% Example code to run RPBM fit on synthetic data
%% Running RPBM noise propagationon simulated 4D data
clc,clear,close all

% Define measurements [ms]
% t = [ 50 100 200 400 ];
t = [ 22, 42, 81, 156, 300 ];
% Compute training data
Ntest = 5e3; diam_factor = 4;
M = Ntest;
bounds = [ 2 2 ; 10 50 ; 0.001 0.05 ]; % D0 in um^2/ms, a in um, kappa in um/ms
[D0_tr,a_tr,kappa_tr,tau_tr,zeta_tr] = RPBMestimation.sample_rpbm_param(bounds,M,diam_factor);
tnorm_all = t./tau_tr;
zeta_all = repmat(zeta_tr,1,length(t));
Dperp_all = RPBMestimation.get_Dt_RPBM_vectorized(tnorm_all, zeta_all);

SNR = 15; Dperp_ref = 0.8; D0_fixed = 2; mask = [];
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
