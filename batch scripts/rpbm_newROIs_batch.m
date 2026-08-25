
clear all
close all

addpath('./supporting_functions_rpbm/');

tic
parent_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell_rescan';
maps_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell_rescan/maps_smoothed';
output_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell_rescan/rpbm_results';
% output_folder = '/Users/gabriellebaxter/Documents/code/rpbm_mapping/d0_greater100ms_slope_outliers';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% roi_folder = '/Users/gabriellebaxter/Documents/segmentation/Dataset001_HEALv1/roi_resampled';
% subjects = dir(roi_folder);
% subjects = {subjects.name};
% subjects = subjects(contains(subjects, 'WCMyofascial2'));

dti_results = readtable(fullfile(parent_folder,'dti_results.csv'));

times = [22,42,81,156,300];
nT = length(times);
nRows = height(dti_results);

% Iterate through rows
for r = 1:nRows
    
    DR = zeros(1, nT); 
    DL = zeros(1, nT); 
    uDR = zeros(1, nT); 
    uDL = zeros(1, nT); 
    for t = 1:nT
        colName = sprintf('RD_median_%04dms', times(t));
        DR(t) = dti_results{r, colName};
        colName = sprintf('AD_median_%04dms', times(t));
        DL(t) = dti_results{r, colName};
        colName = sprintf('RD_std_%04dms', times(t));
        uDR(t) = dti_results{r, colName};
        colName = sprintf('AD_std_%04dms', times(t));
        uDL(t) = dti_results{r, colName};
    end

    UseWeights=1;
    FIX_D0=1;
    if UseWeights==1; W=1./(uDR.^2); else; W=[1,1,1,1,1]; end  

    d0_thresh = 100; % choose how many of the time points to use for Dfix
    Dfix = mean(DL(times>=d0_thresh));                                                 %Fixing D0 to L1(t~inf)
    uDfix=mean(uDL(times>=d0_thresh));

    opt = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[0,0],...
        'Upper',[1000,10],...
        'Startpoint',[100,2],...
        'Robust','bisquare',...
        'Weights',W);
    fitDapp = fittype('D*get_Dt_RPBM(x/tau,zeta)','problem','D','options',opt);
    [fitresult_F,gof_F] = fit(times',DR',fitDapp,'problem',double(Dfix));
    [RPBM_F,uRPBM_F]=RPBM_Process([Dfix,coeffvalues(fitresult_F)],[uDfix,uDfix;confint(fitresult_F)']);%Calculate Various RPBM Parameters and save in structure
    txn=linspace(1,2000,2000);
    RPBM_FIX=squeeze(fitDapp(fitresult_F.tau,fitresult_F.zeta,fitresult_F.D,txn));
    % disp(gof_F)
    
    subjectStr = dti_results.Subject{r};
    roiStr = dti_results.ROI{r};
    disp(subjectStr)
    disp(roiStr)
    % get last 4 characters of Subject
    subjectID = subjectStr(end-3:end);
    
    % build filename
    csv_filename = fullfile(output_folder,sprintf('rpbm_%s_%s.csv', subjectStr, roiStr));
    txt_filename = fullfile(output_folder,sprintf('rpbm_%s_%s.txt', subjectStr, roiStr));

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
toc