
clear all
close all

addpath('./supporting_functions_rpbm/');

parent_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/';

subfolders = dir(fullfile(parent_folder,"MMPS ROI/"));
subjects = {subfolders.name};
subjects = subjects(contains(subjects, '2'));

%for k = 1:length(subjects)
    % vol = subjects{k};
    vol = '2022';
    matching_codes = readtable(fullfile(parent_folder,'matching_codes.csv'));
    % subject_folder = fullfile(parent_folder, ['WCMyofascial' num2str(vol)]);
    subject_folder = fullfile(parent_folder,'C2_013');
    save_folder = fullfile(subject_folder,'rpbm_results_seg2rescan');
    if ~exist(save_folder,"dir") % make rpbm folder if it doesn't exist
        mkdir(save_folder)
    end
    
    % Load ROIs
    % roi_folder = fullfile(parent_folder,'MMPS ROI',vol);
    % roi_file = fullfile(roi_folder,[vol '_mprage.nii']);
    roi_folder = fullfile(subject_folder,'rois');
    roi_file = fullfile(roi_folder,'mprage_reg.nii.gz');
    roi_vol = niftiread(roi_file);
    rois = ["right_masseter","left_masseter","right_temporalis","left_temporalis"];
    
    times = [22,42,81,156,300];
    
    % Get RD and AD values from maps
    results_table = extract_map_values(times,subject_folder,roi_vol,rois);
    writetable(results_table, fullfile(subject_folder,'dti_results.csv'));
    
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
    
        d0_thresh = 0; % choose how many of the time points to use for Dfix
        Dfix = mean(DL(times>=d0_thresh));                                                 %Fixing D0 to L1(t~inf)
        uDfix=uDL(times>=d0_thresh);
        
        opt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0],...
            'Upper',[50,10],...
            'Startpoint',[100,2],...
            'Robust','bisquare',...
            'Weights',W);
        fitDapp = fittype('D*get_Dt_RPBM(x/tau,zeta)','problem','D','options',opt);
        [fitresult_F,gof_F] = fit(times',DR,fitDapp,'problem',double(Dfix));
        [RPBM_F,uRPBM_F]=RPBM_Process([Dfix,coeffvalues(fitresult_F)],[uDfix,uDfix;confint(fitresult_F)']);%Calculate Various RPBM Parameters and save in structure
        txn=linspace(1,2000,2000);
        RPBM_FIX=squeeze(fitDapp(fitresult_F.tau,fitresult_F.zeta,fitresult_F.D,txn));
        disp(gof_F)
    
        csv_filename = fullfile(save_folder, strcat('rpbm_vol', vol, '_', temp_roi,'.csv'));
        txt_filename = fullfile(save_folder, strcat('rpbm_vol', vol, '_', temp_roi, '.txt'));
    
        results = table(txn',RPBM_FIX', 'VariableNames', {'t', 'rpbm_fix'});
        writetable(results, csv_filename);
        
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
%end

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
        if t == 1
            matching_map = '0022ms_AD.nii.gz';
        else
            matching_map = matching_map{end};
        end
        map_array = niftiread(fullfile(maps_folder,matching_map));
        AD_mean = zeros(1,4); AD_std = zeros(1,4); AD_median = zeros(1,4);
        for lbl = 1:4
            AD_mean(lbl) = mean(map_array(roi_vol == lbl), 'omitnan');
            AD_std(lbl) = std(map_array(roi_vol == lbl), 'omitnan');
            AD_median(lbl) = median(map_array(roi_vol == lbl), 'omitnan');
        end
        % RD
        matching_map = matching_map_files(contains(matching_map_files,'RD'));
        if t == 1
            matching_map = '0022ms_RD.nii.gz';
        else
            matching_map = matching_map{end};
        end
        disp(matching_map)
        map_array = niftiread(fullfile(maps_folder,matching_map));
        RD_mean = zeros(1,4); RD_std = zeros(1,4); RD_median = zeros(1,4);
        for lbl = 1:4
            RD_mean(lbl) = mean(map_array(roi_vol == lbl), 'omitnan');
            RD_std(lbl) = std(map_array(roi_vol == lbl), 'omitnan');
            RD_median(lbl) = median(map_array(roi_vol == lbl), 'omitnan');
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
            'VariableNames', ...
            {'diffusion_time','roi', ...
             'RD_mean','RD_std','RD_median','AD_mean','AD_std','AD_median'});
    
        results_table = [results_table; T_new];
    
    
    end
end


