
clear all
close all

addpath('./supporting_functions_rpbm/');

parent_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/';

roi_folder = '/Users/gabriellebaxter/Documents/segmentation/Dataset001_HEALv1/roi_resampled';
subjects = dir(roi_folder);
subjects = {subjects.name};
subjects = subjects(contains(subjects, 'WCMyofascial2'));

%%
for k = 1:2%length(subjects)
    roi_filename = subjects{k};
    vol = roi_filename(13:end-7);
    disp(vol)
    subject_folder = fullfile(parent_folder, ['WCMyofascial' num2str(vol)]);

    save_folder = fullfile(subject_folder,'rpbm_results_negativeslope');
    if ~exist(save_folder,"dir") % make rpbm folder if it doesn't exist
        mkdir(save_folder)
    end

    % Load ROIs
    roi_file = fullfile(roi_folder,roi_filename);
    roi_vol = niftiread(roi_file);
    rois = ["right_masseter","left_masseter","right_temporalis","left_temporalis"];

    times = [22,42,81,156,300];

    % Get RD and AD values from maps
    % results_table = extract_map_values(times,subject_folder,roi_vol,rois);
    % writetable(results_table, fullfile(subject_folder,'dti_results.csv'));

    % Get RD and AD values from maps
    results_table = extract_map_values_negativeslope(times,subject_folder,roi_vol,rois);
    writetable(results_table, fullfile(subject_folder,'dti_results_negativeslope.csv'));

    % Iterate through ROIs and do RPBM
    for i = 1:numel(rois)
        temp_roi = rois(i);

        rows = results_table.roi == temp_roi;
        DL = results_table.AD_median(rows);
        DR = results_table.RD_median(rows);
        uDL = results_table.AD_std(rows);
        uDR = results_table.RD_std(rows);

        UseWeights=1;
        FIX_D0=1;
        if UseWeights==1; W=1./(uDR.^2); else; W=1; end  

        d0_thresh = 200; % choose how many of the time points to use for Dfix
        Dfix = mean(DL(times>=d0_thresh));                                                 %Fixing D0 to L1(t~inf)
        uDfix=uDL(times>=d0_thresh);

        opt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0],...
            'Upper',[1000,10],...
            'Startpoint',[100,2],...
            'Robust','bisquare',...
            'Weights',W);
        fitDapp = fittype('D*get_Dt_RPBM(x/tau,zeta)','problem','D','options',opt);
        [fitresult_F,gof_F] = fit(times',DR,fitDapp,'problem',double(Dfix));
        [RPBM_F,uRPBM_F]=RPBM_Process([Dfix,coeffvalues(fitresult_F)],[uDfix,uDfix;confint(fitresult_F)']);%Calculate Various RPBM Parameters and save in structure
        txn=linspace(1,2000,2000);
        RPBM_FIX=squeeze(fitDapp(fitresult_F.tau,fitresult_F.zeta,fitresult_F.D,txn));
        disp(gof_F)

        % csv_filename = fullfile(save_folder, strcat('rpbm_vol', vol, '_', temp_roi,'.csv'));
        txt_filename = fullfile(save_folder, strcat('rpbm_vol', vol, '_', temp_roi, '.txt'));

        results = table(txn',RPBM_FIX', 'VariableNames', {'t', 'rpbm_fix'});
        % writetable(results, csv_filename);

        % Write params to file
        fileID = fopen(txt_filename, 'w');

        F = fieldnames(RPBM_F);
        fprintf(fileID,'-----  Fitted Results Fixed  ----- \n');

        for i=1:length(F)
            fprintf(fileID,'\n %s = %3.3f +/- %3.3f\n',F{i},RPBM_F.(F{i}),uRPBM_F.(F{i}));
        end
        F = fieldnames(gof_F);
        for i = 1:numel(F)
            fieldName = F{i};
            fieldValue = gof_F.(fieldName);
            fprintf(fileID, '\n %s: %1.4f \n', fieldName, fieldValue);
        end

    end
end

function results_table = extract_map_values(times,subject_folder,roi_vol,rois)
    maps_folder = fullfile(subject_folder,'maps');
    results_table = table();
    for t = 1:length(times)
        dt = times(t);
        disp(dt)
        files = dir(maps_folder);
        names = {files.name};
        matching_map_files = names(contains(names, string(dt)));
        % AD
        matching_map = matching_map_files(contains(matching_map_files,'AD'));
        matching_map = matching_map{end};
        map_array = niftiread(fullfile(maps_folder,matching_map));
        AD_mean = zeros(1,4); AD_std = zeros(1,4); AD_median = zeros(1,4);
        for lbl = 1:4
            AD_mean(lbl) = mean(map_array(roi_vol == lbl), 'omitnan');
            AD_std(lbl) = std(map_array(roi_vol == lbl), 'omitnan');
            AD_median(lbl) = median(map_array(roi_vol == lbl), 'omitnan');
        end
        % RD
        matching_map = matching_map_files(contains(matching_map_files,'RD'));
        matching_map = matching_map{end};
        map_array = niftiread(fullfile(maps_folder,matching_map));
        RD_mean = zeros(1,4); RD_std = zeros(1,4); RD_median = zeros(1,4);
        for lbl = 1:4
            RD_mean(lbl) = mean(map_array(roi_vol == lbl), 'omitnan');
            RD_std(lbl) = std(map_array(roi_vol == lbl), 'omitnan');
            RD_median(lbl) = median(map_array(roi_vol == lbl), 'omitnan');
        end
        % FA
        matching_map = matching_map_files(contains(matching_map_files,'FA'));
        matching_map = matching_map{end};
        map_array = niftiread(fullfile(maps_folder,matching_map));
        map_array = double(map_array);
        FA_mean = zeros(1,4); FA_std = zeros(1,4); FA_median = zeros(1,4);
        for lbl = 1:4
            FA_mean(lbl) = mean(map_array(roi_vol == lbl), 'omitnan');
            FA_std(lbl) = std(map_array(roi_vol == lbl), 'omitnan');
            FA_median(lbl) = median(map_array(roi_vol == lbl), 'omitnan');
        end
    
        % roi = rois(:);   
        roi = categorical(rois(:));
        diffusion_time = repmat(dt, 4, 1);
        
        T_new = table( ...
            diffusion_time, ...
            roi, ...
            RD_mean(:), ...
            RD_std(:), ...
            RD_median(:), ...
            AD_mean(:), ...
            AD_std(:), ...
            AD_median(:), ...
            FA_mean(:), ...
            FA_std(:), ...
            FA_median(:), ...
            'VariableNames', ...
            {'diffusion_time','roi', ...
             'RD_mean','RD_std','RD_median','AD_mean','AD_std','AD_median','FA_mean','FA_std','FA_median'});
    
        results_table = [results_table; T_new];
    
    
    end
end

function results_table = extract_map_values_negativeslope(times,subject_folder,roi_vol,rois)
    maps_folder = fullfile(subject_folder,'maps');
    results_table = table();

    % Load slope map
    slope_array = niftiread(fullfile(subject_folder, 'slopes/RDslope.nii.gz'));
    info = niftiinfo(fullfile(subject_folder, 'slopes/RDr2.nii.gz'));
    info.Datatype = 'uint8';
    slope_array = double(slope_array);

    % Negative slope mask
    neg_mask = slope_array < 0;
    niftiwrite(uint8(neg_mask),fullfile(subject_folder, 'slopes/RD_mask.nii.gz'),info);

    for t = 1:length(times)
        dt = times(t);
        disp(dt)
        files = dir(maps_folder);
        names = {files.name};
        matching_map_files = names(contains(names, string(dt)));
        % AD
        matching_map = matching_map_files(contains(matching_map_files,'AD'));
        matching_map = matching_map{end};
        map_array = niftiread(fullfile(maps_folder,matching_map));
        AD_mean = zeros(1,4); AD_std = zeros(1,4); AD_median = zeros(1,4);
        for lbl = 1:4
            mask = (roi_vol == lbl) & neg_mask;
            AD_mean(lbl)   = mean(map_array(mask), 'omitnan');
            AD_std(lbl)    = std(map_array(mask), 'omitnan');
            AD_median(lbl) = median(map_array(mask), 'omitnan');
        end
        % RD
        matching_map = matching_map_files(contains(matching_map_files,'RD'));
        matching_map = matching_map{end};
        disp(matching_map)
        map_array = niftiread(fullfile(maps_folder,matching_map));
        RD_mean = zeros(1,4); RD_std = zeros(1,4); RD_median = zeros(1,4);
        for lbl = 1:4
            mask = (roi_vol == lbl) & neg_mask;
            RD_mean(lbl)   = mean(map_array(mask), 'omitnan');
            RD_std(lbl)    = std(map_array(mask), 'omitnan');
            RD_median(lbl) = median(map_array(mask), 'omitnan');
        end
        % FA
        matching_map = matching_map_files(contains(matching_map_files,'FA'));
        matching_map = matching_map{end};
        map_array = niftiread(fullfile(maps_folder,matching_map));
        map_array = double(map_array);
        FA_mean = zeros(1,4); FA_std = zeros(1,4); FA_median = zeros(1,4);
        for lbl = 1:4
            mask = (roi_vol == lbl) & neg_mask;
            FA_mean(lbl)   = mean(map_array(mask), 'omitnan');
            FA_std(lbl)    = std(map_array(mask), 'omitnan');
            FA_median(lbl) = median(map_array(mask), 'omitnan');
        end
    
        % roi = rois(:);   
        roi = categorical(rois(:));
        diffusion_time = repmat(dt, 4, 1);
        
        T_new = table( ...
            diffusion_time, ...
            roi, ...
            RD_mean(:), ...
            RD_std(:), ...
            RD_median(:), ...
            AD_mean(:), ...
            AD_std(:), ...
            AD_median(:), ...
            FA_mean(:), ...
            FA_std(:), ...
            FA_median(:), ...
            'VariableNames', ...
            {'diffusion_time','roi', ...
             'RD_mean','RD_std','RD_median','AD_mean','AD_std','AD_median','FA_mean','FA_std','FA_median'});
    
        results_table = [results_table; T_new];
    
    
    end
end

function results_table = extract_map_values_newreg(times,subject_folder,roi_vol,rois)

maps_folder = fullfile(subject_folder,'maps');
slopes_folder = fullfile(subject_folder,'slopes_newreg');

results_table = table();

for t = 1:length(times)

    dt = times(t);
    disp(dt)

    dt_str = sprintf('%04dms',dt);

    AD_mean = zeros(1,4); AD_std = zeros(1,4); AD_median = zeros(1,4);
    RD_mean = zeros(1,4); RD_std = zeros(1,4); RD_median = zeros(1,4);
    FA_mean = zeros(1,4); FA_std = zeros(1,4); FA_median = zeros(1,4);

    for lbl = 1:4

        roi_name = rois{lbl};

        % --- ROI mask ---
        roi_mask = (roi_vol == lbl);

        % --- load slope map for this ROI ---
        slope_file = sprintf('RDslope_%s.nii.gz',roi_name);
        r2_file    = sprintf('RDr2_%s.nii.gz',roi_name);

        slope_array = double(niftiread(fullfile(slopes_folder,slope_file)));
        info = niftiinfo(fullfile(slopes_folder,r2_file));
        info.Datatype = 'uint8';

        % negative slope mask
        neg_mask = slope_array < 0;
        mask = roi_mask & neg_mask;

        % save slope mask
        mask_name = sprintf('RD_mask_%s.nii.gz',roi_name);
        niftiwrite(uint8(mask),fullfile(slopes_folder,mask_name),info);

        % --- construct filenames ---
        if dt == 22
            AD_file = sprintf('%s_AD.nii.gz',dt_str);
            RD_file = sprintf('%s_RD.nii.gz',dt_str);
            FA_file = sprintf('%s_FA.nii.gz',dt_str);
        else
            AD_file = sprintf('%s_AD_%s.nii.gz',dt_str,roi_name);
            RD_file = sprintf('%s_RD_%s.nii.gz',dt_str,roi_name);
            FA_file = sprintf('%s_FA_%s.nii.gz',dt_str,roi_name);
        end

        % --- load maps ---
        AD_array = niftiread(fullfile(maps_folder,AD_file));
        RD_array = niftiread(fullfile(maps_folder,RD_file));
        FA_array = double(niftiread(fullfile(maps_folder,FA_file)));

        % --- statistics ---
        AD_mean(lbl)   = mean(AD_array(mask),'omitnan');
        AD_std(lbl)    = std(AD_array(mask),'omitnan');
        AD_median(lbl) = median(AD_array(mask),'omitnan');

        RD_mean(lbl)   = mean(RD_array(mask),'omitnan');
        RD_std(lbl)    = std(RD_array(mask),'omitnan');
        RD_median(lbl) = median(RD_array(mask),'omitnan');

        FA_mean(lbl)   = mean(FA_array(mask),'omitnan');
        FA_std(lbl)    = std(FA_array(mask),'omitnan');
        FA_median(lbl) = median(FA_array(mask),'omitnan');

    end

    roi = categorical(rois(:));
    diffusion_time = repmat(dt,4,1);

    T_new = table( ...
        diffusion_time,...
        roi,...
        RD_mean(:),...
        RD_std(:),...
        RD_median(:),...
        AD_mean(:),...
        AD_std(:),...
        AD_median(:),...
        FA_mean(:),...
        FA_std(:),...
        FA_median(:),...
        'VariableNames',...
        {'diffusion_time','roi',...
        'RD_mean','RD_std','RD_median',...
        'AD_mean','AD_std','AD_median',...
        'FA_mean','FA_std','FA_median'});

    results_table = [results_table; T_new];

end

end

