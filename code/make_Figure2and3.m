%% Stand-alone: MFCC Figures with MATCHED color scales (TD vs DLD)
% Fixes: ensures TD and DLD heatmaps use the SAME color scale (caxis)
% Affects:
%   Fig1A/Fig1B  (Global MFCC mean|dz| heatmaps)
%   Fig3A/Fig3B  (Representative MFCC dz maps: MFCC-band x coefficient)
%
% Notes:
% - Fig1 uses mean|dz|, so scale is nonnegative.
% - Fig3 uses signed dz, so scale is symmetric around 0.

clear; clc;close all

%% ===== USER SETTINGS =====
dataFolder = pwd;

fileList = { ...
    'FWP_Frnt_stacked.xlsx','FWP_Occ_stacked.xlsx','FWP_Par_stacked.xlsx','FWP_Tmp_stacked.xlsx', ...
    'UWE_Frnt_stacked.xlsx','UWE_Occ_stacked.xlsx','UWE_Par_stacked.xlsx','UWE_Tmp_stacked.xlsx', ...
    'UWR_Frnt_stacked.xlsx','UWR_Occ_stacked.xlsx','UWR_Par_stacked.xlsx','UWR_Tmp_stacked.xlsx'};

taskOrder   = {'FWP','UWE','UWR'};
regionOrder = {'Frnt','Occ','Par','Tmp'};

bandNames = {'delta','theta','alpha','beta','broadband'};
bandPowerCols = [9 10 11 12 13];  %#ok<NASGU> (not used here, but kept for consistency)

% Sheet names (EDIT if needed)
TD_sheetA  = "TD_TherapyA";
TD_sheetB  = "TD_TherapyB";
DLD_sheetA = "DLD_TherapyA";
DLD_sheetB = "DLD_TherapyB";

% Representative conditions for Figure 3 (edit)
repConds = { 'FWP_Frnt', 'UWE_Par', 'FWP_Tmp', 'UWR_Frnt' };

% Save figures?
saveFigs = true;
figOutDir = fullfile(dataFolder, "MFCC_Figures_MatchedColorScale");
if saveFigs && ~exist(figOutDir,'dir'), mkdir(figOutDir); end
%% =========================

%% Build Task–Region keys and labels
rowKeys = cell(numel(taskOrder)*numel(regionOrder),1);
rowLabs = cell(size(rowKeys));
k = 0;
for iT = 1:numel(taskOrder)
    for iR = 1:numel(regionOrder)
        k = k + 1;
        rowKeys{k} = sprintf('%s_%s', taskOrder{iT}, regionOrder{iR});
        rowLabs{k} = sprintf('%s–%s', taskOrder{iT}, regionOrder{iR});
    end
end

%% File lookup map
lut = containers.Map;
for i = 1:numel(fileList)
    f = fileList{i};
    tok = regexp(f, '^(FWP|UWE|UWR)_(Frnt|Occ|Par|Tmp)_stacked\.xlsx$', 'tokens', 'once');
    if ~isempty(tok)
        key = sprintf('%s_%s', tok{1}, tok{2});
        lut(key) = fullfile(dataFolder, f);
    end
end

%% MFCC column blocks (5 MFCC-bands × 13 coeffs)
mfccBlocks = struct();
mfccBlocks.delta     = 14:26;
mfccBlocks.theta     = 27:39;
mfccBlocks.alpha     = 40:52;
mfccBlocks.beta      = 53:65;
mfccBlocks.broadband = 66:78;
mfccBandKeys = {'delta','theta','alpha','beta','broadband'};

%% Compute MFCC dz matrices for TD and DLD
MFCC_dz_TD  = compute_mfcc_dz(lut, rowKeys, mfccBlocks, mfccBandKeys, TD_sheetA,  TD_sheetB);
MFCC_dz_DLD = compute_mfcc_dz(lut, rowKeys, mfccBlocks, mfccBandKeys, DLD_sheetA, DLD_sheetB);

%% ===================== FIGURE 1 (Option 1B) =====================
% Global MFCC therapy sensitivity heatmap: mean |dz| across MFCC-bands
global_TD  = squeeze(mean(abs(MFCC_dz_TD),  2, 'omitnan'));   % (12 x 13)
global_DLD = squeeze(mean(abs(MFCC_dz_DLD), 2, 'omitnan'));   % (12 x 13)

% ---- MATCHED color scale across TD and DLD for Fig1 ----
allGlobal = [global_TD(:); global_DLD(:)];
gMin = min(allGlobal, [], 'omitnan');
gMax = max(allGlobal, [], 'omitnan');
if isempty(gMin) || isempty(gMax) || isnan(gMin) || isnan(gMax) || gMin==gMax
    gMin = 0; gMax = 1; % fallback
end

figure('Name','Fig1A_TD_Global_MFCC_Heatmap_Matched','Color','w');
imagesc(global_TD);
cb = colorbar; ylabel(cb,'mean |dz|');
caxis([gMin gMax]);
title('TD: Global MFCC Therapy Sensitivity (mean |dz| across MFCC-bands)');
xlabel('MFCC coefficient'); ylabel('Task–Region');
set(gca,'YTick',1:12,'YTickLabel',rowLabs,'XTick',1:13);
axis tight;

figure('Name','Fig1B_DLD_Global_MFCC_Heatmap_Matched','Color','w');
imagesc(global_DLD);
cb = colorbar; ylabel(cb,'mean |dz|');
caxis([gMin gMax]);
title('DLD: Global MFCC Therapy Sensitivity (mean |dz| across MFCC-bands)');
xlabel('MFCC coefficient'); ylabel('Task–Region');
set(gca,'YTick',1:12,'YTickLabel',rowLabs,'XTick',1:13);
axis tight;

%% ===================== FIGURE 2 (Band-specific trends) =====================
% (Not a heatmap, so no shared caxis needed)
figure('Name','Fig2_BandSpecific_MFCC_Trends','Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
for b = 1:5
    nexttile;
    profTD  = squeeze(mean(MFCC_dz_TD(:, b, :), 1, 'omitnan'));   % 1x13
    profDLD = squeeze(mean(MFCC_dz_DLD(:, b, :),1, 'omitnan'));   % 1x13

    plot(1:13, profTD, '-o', 'LineWidth', 1.2); hold on;
    plot(1:13, profDLD,'-o', 'LineWidth', 1.2);
    yline(0,'--');

    title(sprintf('MFCC-%s', bandNames{b}));
    xlabel('MFCC coefficient'); ylabel('Mean dz (B−A)');
    legend({'TD','DLD'},'Location','best');
    grid on;
end
nexttile; axis off;

%% ===================== FIGURE 3 (Representative maps) =====================
[repIdx, repLabels] = find_rep_indices(rowKeys, repConds);

% Build all matrices for shared scale
repMatsTD  = cell(numel(repIdx),1);
repMatsDLD = cell(numel(repIdx),1);
for i = 1:numel(repIdx)
    repMatsTD{i}  = squeeze(MFCC_dz_TD(repIdx(i),  :, :));  % 5x13
    repMatsDLD{i} = squeeze(MFCC_dz_DLD(repIdx(i), :, :));  % 5x13
end

% ---- MATCHED (symmetric) color scale across TD and DLD for Fig3 ----
allRep = [];
for i = 1:numel(repIdx)
    allRep = [allRep; repMatsTD{i}(:); repMatsDLD{i}(:)]; %#ok<AGROW>
end
m = max(abs(allRep), [], 'omitnan');
if isempty(m) || isnan(m) || m==0
    m = 1; % fallback
end
repCLim = [-m m];

% figure('Name','Fig3A_TD_Representative_MFCC_Matched','Color','w');
% tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
% for i = 1:numel(repIdx)
%     nexttile;
%     imagesc(repMatsTD{i});
%     cb = colorbar; ylabel(cb,'dz (B−A)');
%     caxis(repCLim);
%     title(['TD: ' strrep(repLabels{i},'_','–')]);
%     xlabel('MFCC coefficient'); ylabel('MFCC-band');
%     set(gca,'YTick',1:5,'YTickLabel',bandNames,'XTick',1:13);
%     axis tight;
% end
% 
% figure('Name','Fig3B_DLD_Representative_MFCC_Matched','Color','w');
% tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
% for i = 1:numel(repIdx)
%     nexttile;
%     imagesc(repMatsDLD{i});
%     cb = colorbar; ylabel(cb,'dz (B−A)');
%     caxis(repCLim);
%     title(['DLD: ' strrep(repLabels{i},'_','–')]);
%     xlabel('MFCC coefficient'); ylabel('MFCC-band');
%     set(gca,'YTick',1:5,'YTickLabel',bandNames,'XTick',1:13);
%     axis tight;
% end

%% Save figures
if saveFigs
    saveas(findobj('Name','Fig1A_TD_Global_MFCC_Heatmap_Matched'), fullfile(figOutDir,'Fig_MFCC1A_TD_GlobalHeatmap_Matched.png'));
    saveas(findobj('Name','Fig1B_DLD_Global_MFCC_Heatmap_Matched'), fullfile(figOutDir,'Fig_MFCC1B_DLD_GlobalHeatmap_Matched.png'));
    saveas(findobj('Name','Fig2_BandSpecific_MFCC_Trends'),        fullfile(figOutDir,'Fig_MFCC2_BandSpecificTrends.png'));
    saveas(findobj('Name','Fig3A_TD_Representative_MFCC_Matched'), fullfile(figOutDir,'Fig_MFCC3A_TD_Representative_Matched.png'));
    saveas(findobj('Name','Fig3B_DLD_Representative_MFCC_Matched'),fullfile(figOutDir,'Fig_MFCC3B_DLD_Representative_Matched.png'));
    fprintf('Saved matched-scale MFCC figures to: %s\n', figOutDir);
end

%% ===================== FUNCTIONS =====================

function MFCC_dz = compute_mfcc_dz(lut, rowKeys, mfccBlocks, mfccBandKeys, sheetA, sheetB)
% MFCC_dz: 12 x 5 x 13 (rows=task-region, mfccBand, coeff)
    MFCC_dz = nan(numel(rowKeys), 5, 13);

    for i = 1:numel(rowKeys)
        key = rowKeys{i};
        if ~isKey(lut, key) || ~isfile(lut(key))
            warning("Missing file for %s", key);
            continue;
        end
        fname = lut(key);

        A = readmatrix(fname, 'Sheet', sheetA);
        B = readmatrix(fname, 'Sheet', sheetB);

        for mb = 1:5
            cols = mfccBlocks.(mfccBandKeys{mb}); % 1x13
            for c = 1:13
                MFCC_dz(i, mb, c) = dz_from_participant_means(A, B, cols(c));
            end
        end
    end
end

function dz = dz_from_participant_means(A, B, valueCol)
% dz = mean(B-A)/std(B-A) using per-participant means
    [uA, mA] = participant_means(A, valueCol);
    [uB, mB] = participant_means(B, valueCol);
    [u, ia, ib] = intersect(uA, uB);
    if numel(u) < 2
        dz = NaN; return;
    end
    d = mB(ib) - mA(ia);
    sd = std(d, 0, 'omitnan');
    if sd == 0 || isnan(sd)
        dz = NaN;
    else
        dz = mean(d,'omitnan') / sd;
    end
end

function [uIDs, means] = participant_means(M, valueCol)
    ids = M(:,1);
    uIDs = unique(ids);
    means = nan(size(uIDs));
    for k = 1:numel(uIDs)
        means(k) = mean(M(ids==uIDs(k), valueCol), 'omitnan');
    end
end

function [idx, labels] = find_rep_indices(rowKeys, repConds)
    idx = nan(1,numel(repConds));
    labels = repConds;
    for i = 1:numel(repConds)
        j = find(strcmp(rowKeys, repConds{i}), 1);
        if isempty(j)
            error("Representative condition not found in rowKeys: %s", repConds{i});
        end
        idx(i) = j;
    end
end
