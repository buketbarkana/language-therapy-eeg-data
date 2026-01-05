%% Figure 1A/1B: Therapy style effects on band power across task × region
% Output:
%   Figure 1A: TD heatmap (effect size dz for TherapyB - TherapyA)
%   Figure 1B: DLD heatmap (effect size dz for TherapyB - TherapyA)
%
% Requires: 12 files in the current folder (or provide folder path):
%   FWP_Frnt_stacked.xlsx ... UWR_Tmp_stacked.xlsx
%
% KEY UPDATE:
%   Both heatmaps use the SAME color range (shared caxis), symmetric about 0.

clear; clc; close all

%% ===== USER SETTINGS =====
dataFolder = pwd;  % change if needed
fileList = { ...
    'FWP_Frnt_stacked.xlsx','FWP_Occ_stacked.xlsx','FWP_Par_stacked.xlsx','FWP_Tmp_stacked.xlsx', ...
    'UWE_Frnt_stacked.xlsx','UWE_Occ_stacked.xlsx','UWE_Par_stacked.xlsx','UWE_Tmp_stacked.xlsx', ...
    'UWR_Frnt_stacked.xlsx','UWR_Occ_stacked.xlsx','UWR_Par_stacked.xlsx','UWR_Tmp_stacked.xlsx'};

% Band-power columns (your layout)
bandCols  = [9 10 11 12 13];                       % delta theta alpha beta broadband
bandNames = {'delta','theta','alpha','beta','broadband'};

taskOrder   = {'FWP','UWE','UWR'};
regionOrder = {'Frnt','Occ','Par','Tmp'};

%% ===== BUILD CONDITION LIST (rows) =====
condLabels = cell(numel(taskOrder)*numel(regionOrder),1);
k = 1;
for t = 1:numel(taskOrder)
    for r = 1:numel(regionOrder)
        condLabels{k} = sprintf('%s-%s', taskOrder{t}, regionOrder{r});
        k = k + 1;
    end
end

TD_dz  = nan(numel(condLabels), numel(bandNames));
DLD_dz = nan(numel(condLabels), numel(bandNames));

%% ===== Helper functions =====
% Find a sheet name ignoring spaces/case
findSheet = @(sheets, target) sheets{ find(strcmp( ...
    upper(regexprep(sheets,'\s+','')), upper(regexprep(target,'\s+',''))), 1) };

% Paired dz for within-subject differences
pairedDz = @(d) mean(d,'omitnan') ./ (std(d,0,'omitnan') + eps);

%% ===== MAIN LOOP OVER FILES =====
for iFile = 1:numel(fileList)

    fn = fileList{iFile};
    fp = fullfile(dataFolder, fn);

    tok = regexp(fn, '^(FWP|UWE|UWR)_(Frnt|Occ|Par|Tmp)_stacked\.xlsx$', 'tokens', 'once');
    if isempty(tok)
        warning('Skipping file (name pattern mismatch): %s', fn);
        continue
    end
    task   = tok{1};
    region = tok{2};

    rowName = sprintf('%s-%s', task, region);
    rowIdx  = find(strcmp(condLabels, rowName), 1);

    [~, sheets] = xlsfinfo(fp);

    %% ===== TD sheets =====
    try
        sh_TD_A = findSheet(sheets, 'TD_TherapyA');
        sh_TD_B = findSheet(sheets, 'TD_TherapyB');

        A = readtable(fp, 'Sheet', sh_TD_A, 'VariableNamingRule','preserve');
        B = readtable(fp, 'Sheet', sh_TD_B, 'VariableNamingRule','preserve');

        pidA = A{:,1}; pidB = B{:,1};
        common = intersect(unique(pidA), unique(pidB));

        for b = 1:numel(bandCols)
            c = bandCols(b);

            xA = nan(numel(common),1);
            xB = nan(numel(common),1);

            for p = 1:numel(common)
                pid = common(p);
                xA(p) = mean(A{pidA==pid, c}, 'omitnan');
                xB(p) = mean(B{pidB==pid, c}, 'omitnan');
            end

            d = xB - xA;                 % TherapyB - TherapyA
            TD_dz(rowIdx, b) = pairedDz(d);
        end
    catch ME
        warning('TD read failed for %s: %s', fn, ME.message);
    end

    %% ===== DLD sheets =====
    try
        sh_DLD_A = findSheet(sheets, 'DLD_TherapyA');
        sh_DLD_B = findSheet(sheets, 'DLD_TherapyB');

        A = readtable(fp, 'Sheet', sh_DLD_A, 'VariableNamingRule','preserve');
        B = readtable(fp, 'Sheet', sh_DLD_B, 'VariableNamingRule','preserve');

        pidA = A{:,1}; pidB = B{:,1};
        common = intersect(unique(pidA), unique(pidB));

        for b = 1:numel(bandCols)
            c = bandCols(b);

            xA = nan(numel(common),1);
            xB = nan(numel(common),1);

            for p = 1:numel(common)
                pid = common(p);
                xA(p) = mean(A{pidA==pid, c}, 'omitnan');
                xB(p) = mean(B{pidB==pid, c}, 'omitnan');
            end

            d = xB - xA;                 % TherapyB - TherapyA
            DLD_dz(rowIdx, b) = pairedDz(d);
        end
    catch ME
        warning('DLD read failed for %s: %s', fn, ME.message);
    end
end

%% ===== SHARED COLOR RANGE (same for TD and DLD) =====
allDz = [TD_dz(:); DLD_dz(:)];
mx = max(abs(allDz), [], 'omitnan');
if isempty(mx) || ~isfinite(mx) || mx==0
    mx = 1; % fallback
end
sharedCLim = [-mx mx];
fprintf('Shared color scale (dz): [%.3f, %.3f]\n', sharedCLim(1), sharedCLim(2));

%% ===== PLOT FIGURE 1A (TD) =====
figure('Color','w','Name','Figure 3A - TD');
imagesc(TD_dz);
axis tight;
set(gca,'YDir','normal'); % optional: top row at top
caxis(sharedCLim);        % <<< SAME RANGE
colorbar;
title('Figure 1A. TD: Therapy effect sizes (dz) for TherapyB − TherapyA');
xticks(1:numel(bandNames)); xticklabels(bandNames);
yticks(1:numel(condLabels)); yticklabels(condLabels);
xlabel('EEG band'); ylabel('Task–Region condition');
set(gca,'TickLabelInterpreter','none');

%% ===== PLOT FIGURE 1B (DLD) =====
figure('Color','w','Name','Figure 3B - DLD');
imagesc(DLD_dz);
axis tight;
set(gca,'YDir','normal'); % optional
caxis(sharedCLim);        % <<< SAME RANGE
colorbar;
title('Figure 1B. DLD: Therapy effect sizes (dz) for TherapyB − TherapyA');
xticks(1:numel(bandNames)); xticklabels(bandNames);
yticks(1:numel(condLabels)); yticklabels(condLabels);
xlabel('EEG band'); ylabel('Task–Region condition');
set(gca,'TickLabelInterpreter','none');

%% Optional: save as PNG (300 dpi)
% exportgraphics(findobj('Name','Figure 3A - TD'),  'Figure3A_TD_therapy_heatmap.png',  'Resolution', 300);
% exportgraphics(findobj('Name','Figure 3B - DLD'), 'Figure3B_DLD_therapy_heatmap.png', 'Resolution', 300);
