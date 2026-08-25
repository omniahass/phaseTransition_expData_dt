%% RPBM fits whole ROI
clc,clear,close all

parent_directory = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell/';
nifti_folder = fullfile(parent_directory,'dwi_smoothed');
bvecbval_folder = fullfile(parent_directory,'bvecbval');
maps_folder = fullfile(parent_directory,'maps_smoothed');
sigma_folder = fullfile(parent_directory,'sigma');
slopes_folder = fullfile(parent_directory,'slopes');
mask_folder = '/Users/gabriellebaxter/Documents/segmentation/Dataset001_HEALv1/roi_resampled';
output_folder = fullfile(parent_directory,'rpbm_maps_bayesian_wholeroi');

sigmas_csv = '/Users/gabriellebaxter/Documents/code/june_reorg/sigma_residuals.csv';

fid = fopen(sigmas_csv,'r');

% skip header
fgetl(fid);

% Subject, ROI, sigma
C = textscan(fid,'%s%s%f', ...
    'Delimiter', ',', ...
    'EndOfLine', '\n');

fclose(fid);

sigma_tbl = table( ...
    string(C{1}), ...
    string(C{2}), ...
    C{3}, ...
    'VariableNames', {'Subject','ROI','sigma'});

subjects = dir(mask_folder);
subjects = {subjects.name};
subjects = subjects(contains(subjects, 'WCMyofascial20'));
baseNames = cellfun(@(x) split(x, '_'), subjects, 'UniformOutput', false);
baseNames = cellfun(@(x) x{1}, baseNames, 'UniformOutput', false);
uniqueNames = unique(baseNames);

outFile = fullfile(parent_directory, 'rpbm_results_bayesian_wholeroi_sigmaresidual.csv');
if exist(outFile, 'file')
    delete(outFile);
end
% subjects = subjects(contains(subjects, 'C2_'));
rois = ["left_masseter", "right_masseter", ...
        "left_temporalis", "right_temporalis"];
%%
for i = 1:length(uniqueNames)
    subject = uniqueNames{i};
    subject = subject(1:end-7);
    disp(subject)
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
        if t == 1
            filename = fullfile(maps_folder, ...
                sprintf('%s_%s_RD.nii.gz', subject, times{t}));
        else
            filename = fullfile(maps_folder, ...
                sprintf('%s_%s_RD_withinscanreg.nii.gz', subject, times{t}));
        end
        nii = load_untouch_nii(filename);
        vol = nii.img;
        if t == 1
            [nx, ny, nz] = size(vol);
            Dperp = zeros(nx, ny, nz, nT, 'single');  % preallocate
        end
        Dperp(:,:,:,t) = vol;
    end
    % Load sigma
    filename = fullfile(sigma_folder, sprintf('%s_sigma.nii', subject));
    nii_sigma = load_untouch_nii(filename); 
    info = niftiinfo(filename); 
    % Load slopes
    filename = fullfile(slopes_folder, sprintf('%s_RD_slope.nii.gz', subject));
    nii_slope = load_untouch_nii(filename); 
    slope = nii_slope.img;
    % Define a mask where you want the RPBM fits
    filename = fullfile(mask_folder, sprintf('%s.nii.gz', subject));
    nii_mask = load_untouch_nii(filename); 
    mask = nii_mask.img;
    mask_new = mask;
    mask_new(slope >= 0) = 0;
    
    labels = [1 2 3 4];
    
    t = [22, 42, 81, 156, 300];
    
    sigma_b0_norm = nii_sigma.img ./ (b0 + 1e-6);
    
    Dref = [1.7 0.9];
    mask_for_fit = mask_new > 0;
    [ sigma_Dpar , sigma_Dperp ] = ...
        RPBMestimation.LLS_noiseprop_DTI_numerical(Dref,bvals,dirs,sigma_b0_norm,mask_for_fit);
    
    D0_fixed = 2;
    PRoptions.Ntrain = 2e4;
    PRoptions.diam_factor = 6.29;
    PRoptions.bounds = [ D0_fixed D0_fixed ; 10 60 ; 0.001 0.05 ];
    PRoptions.PRdegree = 4;
    
    results = struct();
    
    for L = labels
        roi_name = rois(L);
        
        % --- Create binary mask for this label ---
        roi_mask = (mask_new == L);
        
        if nnz(roi_mask) == 0
            fprintf('Label %d is empty, skipping\n', L);
            continue;
        end
        
        % --- Average Dperp over ROI ---
        Dperp_roi = zeros(1, length(t));
        
        for ti = 1:length(t)
            vol = Dperp(:,:,:,ti);
            Dperp_roi(ti) = mean(vol(roi_mask), 'omitnan');
        end
        
        % --- Average noise ---
        % N = nnz(roi_mask);
        % sigma_vox = sigma_Dperp(roi_mask);
        % sigma_roi = 0.5*median(sigma_vox,'omitnan');

        % --- Get sigma from CSV ---
        match_idx = sigma_tbl.Subject == string(subject) & ...
            sigma_tbl.ROI == string(roi_name);

        if ~any(match_idx)
            warning('No sigma found for %s / %s', subject, roi_name);
            continue;
        end
        
        sigma_roi = sigma_tbl.sigma(find(match_idx,1));
        
        % --- Reshape for fitting (needs 4D input) ---
        Dperp_fit = reshape(Dperp_roi, [1 1 1 length(t)]);
        sigma_fit = reshape(sigma_roi, [1 1 1]);
        mask_fit = 1;  % single voxel
        
        % --- Run RPBM ---
        % sigma_fit=1;
        % sigma_fit = sigma_fit*10;
        param = RPBMestimation.fit4Ddata( ...
            Dperp_fit, t, sigma_fit, D0_fixed, mask_fit, PRoptions);
        
        % --- Store ---
        results(L).Dperp = Dperp_roi;
        results(L).sigma = sigma_roi;
        results(L).param = param;

        row = table();

        row.Subject = string(subject);
        row.ROI     = roi_name;
        
        row.Dperp_22  = Dperp_roi(1);
        row.Dperp_42  = Dperp_roi(2);
        row.Dperp_81  = Dperp_roi(3);
        row.Dperp_156 = Dperp_roi(4);
        row.Dperp_300 = Dperp_roi(5);
        
        row.sigma_roi = sigma_roi;
        
        row.a     = param.a;
        row.kappa = param.kappa;
        
        writetable(row, outFile, 'WriteMode','append');
                
        fprintf('Finished label %s\n', roi_name);
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
