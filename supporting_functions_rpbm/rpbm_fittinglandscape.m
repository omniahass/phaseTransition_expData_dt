
clear all
close all

addpath('./supporting_functions_rpbm/');

parent_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell';
maps_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell/maps_smoothed';
% output_folder = '/Users/gabriellebaxter/Documents/HEAL_Volunteers/cornell/rpbm_results_negativeslope';
output_folder = '/Users/gabriellebaxter/Documents/code/rpbm_mapping/t100z0p5';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%%

% roi_folder = '/Users/gabriellebaxter/Documents/segmentation/Dataset001_HEALv1/roi_resampled';
% subjects = dir(roi_folder);
% subjects = {subjects.name};
% subjects = subjects(contains(subjects, 'WCMyofascial2'));

dti_results = readtable(fullfile(parent_folder,'dti_results_maskedreg.csv'));

times = [22,42,81,156,300];
nT = length(times);
nRows = height(dti_results);

% Iterate through rows
for r = 1%:nRows
    
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
    if UseWeights==1; W=1./(uDR.^2); else; W=1; end  

    d0_thresh = 200; % choose how many of the time points to use for Dfix
    Dfix = mean(DL(times>=d0_thresh));                                                 %Fixing D0 to L1(t~inf)
    uDfix=uDL(times>=d0_thresh);
    % 
    % opt = fitoptions('Method','NonlinearLeastSquares',...
    %     'Lower',[0,0],...
    %     'Upper',[1000,10],...
    %     'Startpoint',[100,0.5],...
    %     'Robust','bisquare',...
    %     'Weights',W);
    % fitDapp = fittype('D*get_Dt_RPBM(x/tau,zeta)','problem','D','options',opt);
    % [fitresult_F,gof_F] = fit(times',DR',fitDapp,'problem',double(Dfix));
    % [RPBM_F,uRPBM_F]=RPBM_Process([Dfix,coeffvalues(fitresult_F)],[uDfix,uDfix;confint(fitresult_F)']);%Calculate Various RPBM Parameters and save in structure
    
    % --- Define multiple starting points ---
    startPoints = [
        10   0.1
        50   0.5
        100  1
        200  2
        500  5
    ];
    
    nStarts = size(startPoints,1);
    
    tau_all   = nan(nStarts,1);
    zeta_all  = nan(nStarts,1);
    r2_all    = nan(nStarts,1);
    sse_all   = nan(nStarts,1);
    
    bestFit = [];
    bestGof = [];
    bestR2 = -inf;
    
    for s = 1:nStarts
    
        fprintf('Trying startpoint [%g, %g]\n', ...
            startPoints(s,1), startPoints(s,2));
    
        opt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0],...
            'Upper',[1000,10],...
            'Startpoint',startPoints(s,:),...
            'Robust','bisquare',...
            'Weights',W);
    
        fitDapp = fittype( ...
            'D*get_Dt_RPBM(x/tau,zeta)', ...
            'problem','D', ...
            'options',opt);
    
        try
            [fitresult_tmp, gof_tmp] = fit( ...
                times', DR', fitDapp, ...
                'problem', double(Dfix));
    
            % Store fit results
            tau_all(s)  = fitresult_tmp.tau;
            zeta_all(s) = fitresult_tmp.zeta;
            r2_all(s)   = gof_tmp.rsquare;
            sse_all(s)  = gof_tmp.sse;
    
            % Keep best R² fit
            if gof_tmp.rsquare > bestR2
                bestR2 = gof_tmp.rsquare;
                bestFit = fitresult_tmp;
                bestGof = gof_tmp;
            end
    
        catch ME
            warning('Fit failed for startpoint [%g %g]: %s', ...
                startPoints(s,1), ...
                startPoints(s,2), ...
                ME.message);
        end
    end
    
    % Final chosen fit
    fitresult_F = bestFit;
    gof_F = bestGof;
    
    % Display table of results
    disp(table( ...
        startPoints(:,1), ...
        startPoints(:,2), ...
        tau_all, ...
        zeta_all, ...
        r2_all, ...
        sse_all, ...
        'VariableNames', ...
        {'tau_start','zeta_start','tau_fit', ...
         'zeta_fit','R2','SSE'}))
    
    % Variation summary
    fprintf('\nVariation across starts:\n');
    
    fprintf('tau: mean = %.3f, std = %.3f, range = [%.3f %.3f]\n', ...
        mean(tau_all,'omitnan'), ...
        std(tau_all,'omitnan'), ...
        min(tau_all), ...
        max(tau_all));
    
    fprintf('zeta: mean = %.3f, std = %.3f, range = [%.3f %.3f]\n', ...
        mean(zeta_all,'omitnan'), ...
        std(zeta_all,'omitnan'), ...
        min(zeta_all), ...
        max(zeta_all));
    
    fprintf('R²: mean = %.4f, std = %.4f, range = [%.4f %.4f]\n', ...
        mean(r2_all,'omitnan'), ...
        std(r2_all,'omitnan'), ...
        min(r2_all), ...
        max(r2_all));
    
    fprintf('\nBest fit:\n');
    fprintf('tau = %.3f\n', fitresult_F.tau);
    fprintf('zeta = %.3f\n', fitresult_F.zeta);
    fprintf('R² = %.4f\n', gof_F.rsquare);
    fprintf('SSE = %.4f\n', gof_F.sse);
        
    % Process RPBM outputs
    [RPBM_F,uRPBM_F] = RPBM_Process( ...
        [Dfix, coeffvalues(fitresult_F)], ...
        [uDfix, uDfix; confint(fitresult_F)']);
    
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
