%% Stand-alone code: compute dz matrices (TD & DLD) + generate Table 3
% Table 3. Regional therapy sensitivity across tasks and frequency bands
% - Computes dz = mean(B-A)/std(B-A) for each Task×Region×Band cell
% - Builds Table 3 by summarizing within each Region and Group
%
% Assumptions:
%   - 12 stacked files exist in dataFolder (listed in fileList)
%   - Each file has TWO sheets for each group:
%       TD_sheetA / TD_sheetB
%       DLD_sheetA / DLD_sheetB
%   - Columns: Col1=participant, Col2=word_id (ignored)
%   - Band power columns: delta=9, theta=10, alpha=11, beta=12, broadband=13

clear; clc;close all;

%% ===== USER SETTINGS =====
dataFolder = pwd;  % change if needed

fileList = { ...
    'FWP_Frnt_stacked.xlsx','FWP_Occ_stacked.xlsx','FWP_Par_stacked.xlsx','FWP_Tmp_stacked.xlsx', ...
    'UWE_Frnt_stacked.xlsx','UWE_Occ_stacked.xlsx','UWE_Par_stacked.xlsx','UWE_Tmp_stacked.xlsx', ...
    'UWR_Frnt_stacked.xlsx','UWR_Occ_stacked.xlsx','UWR_Par_stacked.xlsx','UWR_Tmp_stacked.xlsx'};

bandCols  = [9 10 11 12 13];
bandNames = {'delta','theta','alpha','beta','broadband'};

taskOrder   = {'FWP','UWE','UWR'};
regionOrder = {'Frnt','Occ','Par','Tmp'};

% Sheet names (EDIT to match your Excel)
TD_sheetA  = "TD_TherapyA";
TD_sheetB  = "TD_TherapyB";
DLD_sheetA = "DLD_TherapyA";
DLD_sheetB = "DLD_TherapyB";

% Table 3 rules
taskRankingMode = "meanAbs";   % "meanAbs" or "maxAbs"
bandRankingMode = "meanAbs";   % "meanAbs" or "maxAbs"
topKTasks = 2;
topKBands = 2;

% Export
exportExcel = true;
outXlsx = "Table3_RegionalTherapySensitivity.xlsx";
%% =========================

%% Build row labels & keys (12 rows: task-region)
rowKeys = cell(numel(taskOrder)*numel(regionOrder),1); % e.g., FWP_Frnt
idx = 0;
for iT = 1:numel(taskOrder)
    for iR = 1:numel(regionOrder)
        idx = idx + 1;
        rowKeys{idx} = sprintf('%s_%s', taskOrder{iT}, regionOrder{iR});
    end
end

%% Compute dz matrices: rows=12 task-region, cols=5 bands
dz_TD  = compute_dz_matrix(dataFolder, fileList, rowKeys, bandCols, TD_sheetA,  TD_sheetB);
dz_DLD = compute_dz_matrix(dataFolder, fileList, rowKeys, bandCols, DLD_sheetA, DLD_sheetB);

fprintf("TD dz range:  min=%.3f  max=%.3f\n", min(dz_TD(:),[],'omitnan'),  max(dz_TD(:),[],'omitnan'));
fprintf("DLD dz range: min=%.3f  max=%.3f\n", min(dz_DLD(:),[],'omitnan'), max(dz_DLD(:),[],'omitnan'));

%% Prepare Table 3
T3 = build_table3(dz_TD, dz_DLD, bandNames, taskRankingMode, bandRankingMode, topKTasks, topKBands);

disp("Table 3. Regional therapy sensitivity across tasks and frequency bands");
disp(T3);

if exportExcel
    writetable(T3, outXlsx, 'FileType','spreadsheet');
    fprintf("Saved: %s\n", outXlsx);
end

%% ===================== FUNCTIONS =====================

function dzMat = compute_dz_matrix(dataFolder, fileList, rowKeys, bandCols, sheetA, sheetB)
% Returns dzMat: (12 x 5) where rows=task-region, cols=bands

    nCond  = numel(rowKeys);
    nBands = numel(bandCols);
    dzMat = nan(nCond, nBands);

    % File lookup from list
    lut = containers.Map;
    for k = 1:numel(fileList)
        f = fileList{k};
        tok = regexp(f, '^(FWP|UWE|UWR)_(Frnt|Occ|Par|Tmp)_stacked\.xlsx$', 'tokens', 'once');
        if ~isempty(tok)
            key = sprintf('%s_%s', tok{1}, tok{2}); % e.g., FWP_Frnt
            lut(key) = fullfile(dataFolder, f);
        end
    end

    for iC = 1:nCond
        key = rowKeys{iC};

        if ~isKey(lut, key)
            warning("Missing in fileList: %s", key);
            continue;
        end
        fname = lut(key);

        if ~isfile(fname)
            warning("Missing file on disk: %s", fname);
            continue;
        end

        % Read A/B
        A = readmatrix(fname, 'Sheet', sheetA);
        B = readmatrix(fname, 'Sheet', sheetB);

        for b = 1:nBands
            [uA, meanA] = participant_means(A, bandCols(b));
            [uB, meanB] = participant_means(B, bandCols(b));

            [u, ia, ib] = intersect(uA, uB);
            if numel(u) < 2
                dzMat(iC, b) = NaN;
                continue;
            end

            d = meanB(ib) - meanA(ia);  % (B - A)
            sd = std(d, 0, 'omitnan');
            if sd == 0 || isnan(sd)
                dzMat(iC, b) = NaN;
            else
                dzMat(iC, b) = mean(d,'omitnan') / sd;
            end
        end
    end
end

function [uIDs, means] = participant_means(M, valueCol)
% Unique participant IDs and per-participant mean across rows
    ids = M(:,1);
    uIDs = unique(ids);
    means = nan(size(uIDs));
    for k = 1:numel(uIDs)
        means(k) = mean(M(ids==uIDs(k), valueCol), 'omitnan');
    end
end

function T = build_table3(dz_TD, dz_DLD, bandNames, taskMode, bandMode, topKTasks, topKBands)
% Builds Table 3 with 8 rows (4 regions × 2 groups)

    % Row indices per region in the fixed 12-row ordering
    idx.Frnt = [1 5 9];
    idx.Occ  = [2 6 10];
    idx.Par  = [3 7 11];
    idx.Tmp  = [4 8 12];

    regionsFull = {'Frontal','Occipital','Parietal','Temporal'};
    regionKey   = {'Frnt','Occ','Par','Tmp'};

    T = table;
    row = 0;

    for r = 1:numel(regionKey)
        rk = regionKey{r};
        ridx = idx.(rk);

        % TD
        row = row + 1;
        [taskStr, bandStr, rangeStr] = summarize_region(dz_TD, ridx, bandNames, taskMode, bandMode, topKTasks, topKBands);
        T.Region{row,1} = regionsFull{r};
        T.Group{row,1} = 'TD';
        T.Tasks_with_strongest_therapy_effects{row,1} = taskStr;
        T.Dominant_frequency_bands{row,1} = bandStr;
        T.Effect_size_range_dz{row,1} = rangeStr;

        % DLD
        row = row + 1;
        [taskStr, bandStr, rangeStr] = summarize_region(dz_DLD, ridx, bandNames, taskMode, bandMode, topKTasks, topKBands);
        T.Region{row,1} = regionsFull{r};
        T.Group{row,1} = 'DLD';
        T.Tasks_with_strongest_therapy_effects{row,1} = taskStr;
        T.Dominant_frequency_bands{row,1} = bandStr;
        T.Effect_size_range_dz{row,1} = rangeStr;
    end
end

function [taskStr, bandStr, rangeStr] = summarize_region(dzMat, ridx, bandNames, taskMode, bandMode, topKTasks, topKBands)
% dzMat: 12x5; ridx: 1x3 (FWP, UWE, UWR rows for that region)

    tasks = {'FWP','UWE','UWR'};     % fixed ordering within a region
    X = dzMat(ridx, :);             % 3x5

    % Rank tasks
    switch lower(taskMode)
        case "meanabs"
            scoreT = mean(abs(X), 2, 'omitnan');
        case "maxabs"
            scoreT = max(abs(X), [], 2, 'omitnan');
        otherwise
            error("Unknown taskRankingMode");
    end
    [~, ordT] = sort(scoreT, 'descend');
    kT = min(topKTasks, numel(tasks));
    taskStr = strjoin(tasks(ordT(1:kT)), ', ');

    % Rank bands
    switch lower(bandMode)
        case "meanabs"
            scoreB = mean(abs(X), 1, 'omitnan');
        case "maxabs"
            scoreB = max(abs(X), [], 1, 'omitnan');
        otherwise
            error("Unknown bandRankingMode");
    end
    [~, ordB] = sort(scoreB, 'descend');
    kB = min(topKBands, numel(bandNames));
    bandStr = strjoin(bandNames(ordB(1:kB)), ', ');

    % Range across all tasks and bands in that region
    meanAbs = mean(abs(X(:)), 'omitnan');
    xmin = min(X(:), [], 'omitnan');
    xmax = max(X(:), [], 'omitnan');
    rangeStr = sprintf('%.2f – %.2f (mean|dz|=%.2f)', xmin, xmax, meanAbs);
end
