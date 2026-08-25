% Plot ROI results from dti_rpbm_results.csv
%
% Figure 1: RPBM a and kappa per ROI. These are bars, not curves - the RPBM fit
%           uses all five diffusion times at once, so there is no time axis.
% Figure 2: AD and RD against diffusion time, all four ROIs on one axes.
%
% Set save_figures to true to write them into <parent_folder>/figures

clear all
close all
clc

% SET PARENT FOLDER HERE, MAKE SURE IT HAS / ON THE END
% run_heal_pipeline.sh sets the HEAL_PARENT environment variable, otherwise the
% hard-coded default below is used
parent_folder = getenv('HEAL_PARENT');
if isempty(parent_folder)
    %parent_folder = '/mnt/fieree01lab/omnia/projects/HEAL/data/WCMyofascial2075';
    parent_folder = '/Volumes/labspace/omnia/projects/HEAL/data/WCMyofascial2075';
end
if ~endsWith(parent_folder,filesep) % this script needs the trailing /
    parent_folder = [parent_folder filesep];
end
fprintf('Parent folder: %s\n',parent_folder);

% ------------------------------- settings ----------------------------------
diffusion_times = [22 42 81 156 300];                 % ms
time_labels = {'0022ms','0042ms','0081ms','0156ms','0300ms'};
rois = {'right_masseter','left_masseter','right_temporalis','left_temporalis'};
roi_labels = {'R masseter','L masseter','R temporalis','L temporalis'};
params = {'AD','RD'};

statistic = 'mean';        % 'mean' or 'median', the csv holds both
show_errorbars = true;     % std across voxels in the ROI, not an error on the mean
save_figures = false;      % write png and fig into <parent_folder>/figures

markers = {'o','s','^','d'};
roi_colours = [0.85 0.33 0.10;    % R masseter
               0.93 0.69 0.13;    % L masseter
               0.00 0.45 0.74;    % R temporalis
               0.30 0.75 0.93];   % L temporalis

% close shades within each parameter so AD and RD read as two groups
ad_colours = [0.55 0.08 0.08; 0.80 0.25 0.10; 0.93 0.49 0.16; 0.98 0.70 0.35];
rd_colours = [0.04 0.22 0.50; 0.10 0.44 0.75; 0.25 0.65 0.90; 0.55 0.82 0.95];
group_colours = {ad_colours, rd_colours};
group_styles  = {'-','--'};

% ------------------------------- read data ---------------------------------
csv_file = fullfile(parent_folder,'dti_rpbm_results_python.csv');
if ~exist(csv_file,'file')
    error('results file not found: %s', csv_file);
end
T = readtable(csv_file,'VariableNamingRule','preserve');
fprintf('Read %d ROI rows from %s\n', height(T), csv_file);

[~, subject] = fileparts(parent_folder(1:end-1));

% ------------------- figure 1: RPBM parameters per ROI ---------------------
rpbm_params = {'a','kappa'};
rpbm_units  = {'\mum','\mum/ms'};

figure('Position',[100 100 1100 450],'Color','w');

for p = 1:length(rpbm_params)
    param = rpbm_params{p};
    subplot(1,length(rpbm_params),p); hold on

    values = nan(1,length(rois));
    spread = nan(1,length(rois));
    for r = 1:length(rois)
        row = strcmp(T.ROI, rois{r});
        if ~any(row)
            warning('ROI %s not in the csv, skipping', rois{r});
            continue
        end
        values(r) = T.(sprintf('%s_%s',param,statistic))(row);
        spread(r) = T.(sprintf('%s_std',param))(row);
        bar(r, values(r), 0.6, 'FaceColor', roi_colours(r,:), 'EdgeColor','none');
    end

    if show_errorbars
        errorbar(1:length(rois), values, spread, 'k', 'LineStyle','none', ...
            'LineWidth',1.2, 'CapSize',8);
    end

    ylabel(sprintf('%s (%s)',param,rpbm_units{p}))
    title(sprintf('RPBM %s - ROI %s',param,statistic))
    xlim([0.4 length(rois)+0.6])
    set(gca,'FontSize',13,'XTick',1:length(rois),'XTickLabel',roi_labels, ...
        'XTickLabelRotation',20,'Box','off')
    grid on

    fprintf('%-6s %-8s %s\n', param, statistic, num2str(values,'%10.4f'));
end

sgtitle(sprintf('%s - RPBM parameters per ROI', strrep(subject,'_','\_')), ...
    'FontSize',15)
save_figure(save_figures, parent_folder, sprintf('roi_rpbm_%s',statistic));

% ------------- figure 2: AD and RD vs diffusion time, one axes -------------
figure('Position',[100 100 850 600],'Color','w'); hold on
legend_entries = cell(1,length(params)*length(rois));
n = 0;

for p = 1:length(params)
    param = params{p};
    for r = 1:length(rois)
        row = strcmp(T.ROI, rois{r});
        if ~any(row)
            continue
        end

        values = nan(1,length(time_labels));
        for t = 1:length(time_labels)
            values(t) = T.(sprintf('%s_%s_%s',param,statistic,time_labels{t}))(row);
        end

        % plain lines here - 8 series with error bars is unreadable
        plot(diffusion_times, values, ...
            'Color', group_colours{p}(r,:), 'LineStyle', group_styles{p}, ...
            'Marker', markers{r}, 'MarkerSize', 7, ...
            'MarkerFaceColor', group_colours{p}(r,:), 'LineWidth', 1.6);

        n = n + 1;
        legend_entries{n} = sprintf('%s %s', param, roi_labels{r});
        fprintf('%-4s %-8s %-18s %s\n', param, statistic, rois{r}, ...
            num2str(values,'%8.3f'));
    end
end
legend_entries = legend_entries(1:n);

xlabel('Diffusion time (ms)')
ylabel('Diffusivity (\mum^2/ms)')
title(sprintf('%s - AD and RD, ROI %s', strrep(subject,'_','\_'), statistic))
xlim([0 320])
ylim([0.5 2.5])
set(gca,'FontSize',13,'XTick',diffusion_times,'Box','off')
grid on
legend(legend_entries,'Location','eastoutside','Box','off','FontSize',11)

save_figure(save_figures, parent_folder, sprintf('roi_AD_RD_combined_%s',statistic));

% ---------------------------------------------------------------------------
function save_figure(do_save, parent_folder, name)
    if ~do_save
        return
    end
    figure_folder = fullfile(parent_folder,'figures');
    if ~exist(figure_folder,'dir')
        mkdir(figure_folder)
    end
    out_file = fullfile(figure_folder,name);
    saveas(gcf,[out_file '.png']);
    savefig(gcf,[out_file '.fig']);
    fprintf('Saved %s.png\n', out_file);
end
