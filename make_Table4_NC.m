%% Stand-alone MATLAB script to generate Table 5
% Table 5. Region-level strongest therapy effect (minimum p across tasks)
%
% For each REGION and GROUP (TD, DLD):
%   - Evaluate all Task × Band combinations within that region (3 tasks × 5 bands = 15 cells)
%   - Compute paired dz, paired t-test p, and within-task-region FDR q (across 5 bands)
%   - Select the single "strongest" therapy-sensitive condition as the one with MINIMUM p
%       Tie-break #1: larger |dz|
%       Tie-break #2: smaller q
%   - Report: task–band (dir), dz, p (q)
%
% IMPORTANT:
% q-values are computed within each Task–Region condition across the 5 bands,
% consistent with your Table 4 definition. The q reported here is the q-value
% for the selected Task–Band cell.
%
% Assumptions:
% - 12 stacked Excel files exist in dataFolder (fileList)
% - Each file has sheets for TD and DLD TherapyA/B (edit sheet names below)
% - Columns: Col1 participant ID; Col2 word_id (ignored)
% - Band power columns: delta=9, theta=10, alpha=11, beta=12, broadband=13
%
% Output:
% - MATLAB table T5 (4 rows: Frontal/Occipital/Parietal/Temporal)
% - Optional Excel export

clear; clc;

%% ===== USER SETTINGS =====
dataFolder = pwd;

fileList = { ...
    'FWP_Frnt_stacked.xlsx','FWP_Occ_stacked.xlsx','FWP_Par_stacked.xlsx','FWP_Tmp_stacked.xlsx', ...
    'UWE_Frnt_stacked.xlsx','UWE_Occ_stacked.xlsx','UWE_Par_stacked.xlsx','UWE_Tmp_stacked.xlsx', ...
    'UWR_Frnt_stacked.xlsx','UWR_Occ_stacked.xlsx','UWR_Par_stacked.xlsx','UWR_Tmp_stacked.xlsx'};

taskOrder   = {'FWP','UWE','UWR'};                 % task labels
regionOrder = {'Frnt','Occ','Par','Tmp'};          % region keys
regionFull  = {'Frontal','Occipital','Parietal','Temporal'};

bandCols  = [9 10 11 12 13];
bandNames = {'delta','theta','alpha','beta','broadband'};

% Sheet names (EDIT if your files use different names)
TD_sheetA  = "TD_TherapyA";
TD_sheetB  = "TD_TherapyB";
DLD_sheetA = "DLD_TherapyA";
DLD_sheetB = "DLD_TherapyB";

% Export
exportExcel = true;
outXlsx = "Table5_RegionMinP_StrongestTherapyEffect.xlsx";
%% =========================

%% Build file lookup
lut = containers.Map;
for i = 1:numel(fileList)
    f = fileList{i};
    tok = regexp(f, '^(FWP|UWE|UWR)_(Frnt|Occ|Par|Tmp)_stacked\.xlsx$', 'tokens', 'once');
    if ~isempty(tok)
        key = sprintf('%s_%s', tok{1}, tok{2});
        lut(key) = fullfile(dataFolder, f);
    else
        warning("Ignoring unexpected filename: %s", f);
    end
end

%% For each region: find min-p cell across tasks & bands for TD and DLD
T5 = table;
T5.Region = string(regionFull(:));

TD_cell  = strings(4,1);
DLD_cell = strings(4,1);

for r = 1:4
    regKey = regionOrder{r};

    % Collect all cells for TD and DLD across tasks and bands
    cellsTD  = [];
    cellsDLD = [];

    for t = 1:numel(taskOrder)
        task = taskOrder{t};
        key = sprintf('%s_%s', task, regKey);

        if ~isKey(lut, key) || ~isfile(lut(key))
            warning("Missing file for %s", key);
            continue;
        end
        fname = lut(key);

        % Compute dz/p/q vectors (1x5) for this task-region
        [dzTD, pTD, qTD, dirTD, NTD]   = compute_dz_p_q(fname, TD_sheetA,  TD_sheetB,  bandCols);
        [dzDLD,pDLD,qDLD,dirDLD,NDLD] = compute_dz_p_q(fname, DLD_sheetA, DLD_sheetB, bandCols);

        % Add 5 band-cells to pool
        for b = 1:numel(bandNames)
            cellsTD  = [cellsTD;  make_cell(task, regKey, bandNames{b}, dirTD(b),  dzTD(b),  pTD(b),  qTD(b),  NTD)];
            cellsDLD = [cellsDLD; make_cell(task, regKey, bandNames{b}, dirDLD(b), dzDLD(b), pDLD(b), qDLD(b), NDLD)];
        end
    end

    % Pick best (min p) cell for each group in this region
    bestTD  = pick_best_min_p(cellsTD);
    bestDLD = pick_best_min_p(cellsDLD);

    TD_cell(r)  = format_region_best(bestTD);
    DLD_cell(r) = format_region_best(bestDLD);
end

T5.TD  = TD_cell;
T5.DLD = DLD_cell;

disp("Table 5. Region-level strongest therapy effect (minimum p across tasks)");
disp(T5);

if exportExcel
    writetable(T5, outXlsx, 'FileType','spreadsheet');
    fprintf("Saved: %s\n", outXlsx);
end

%% ===================== FUNCTIONS =====================

function row = make_cell(task, region, band, arrow, dz, p, q, N)
% Store one task-region-band cell as a struct row
    row = struct();
    row.task = string(task);
    row.region = string(region);
    row.band = string(band);
    row.arrow = string(arrow);
    row.dz = dz;
    row.p = p;
    row.q = q;
    row.N = N;
end

function best = pick_best_min_p(cells)
% Select cell with minimum p; tie-break by larger |dz| then smaller q
    if isempty(cells)
        best = make_cell("NA","NA","NA","?",NaN,NaN,NaN,NaN);
        return;
    end

    pvals = arrayfun(@(c) c.p, cells);
    dzvals = arrayfun(@(c) c.dz, cells);
    qvals = arrayfun(@(c) c.q, cells);

    valid = ~isnan(pvals);
    if ~any(valid)
        best = cells(1);
        return;
    end

    pvals(~valid) = Inf; % ignore NaN p
    minP = min(pvals);

    cand = find(pvals == minP);
    if numel(cand) == 1
        best = cells(cand);
        return;
    end

    % tie-break: max |dz|
    absdz = abs(dzvals(cand));
    maxAbs = max(absdz);
    cand2 = cand(absdz == maxAbs);

    if numel(cand2) == 1
        best = cells(cand2);
        return;
    end

    % tie-break: min q
    qsub = qvals(cand2);
    qsub(isnan(qsub)) = Inf;
    [~, j] = min(qsub);
    best = cells(cand2(j));
end

function s = format_region_best(best)
% Output style: "FWP–alpha (↑), dz=0.706, p=0.086 (q=0.199) [N=8]"
    if isnan(best.dz) || isnan(best.p) || isnan(best.q)
        s = "NA";
        return;
    end
    s = sprintf('%s–%s (%s), dz=%.3f, p=%.3f (q=%.3f) [N=%d]', ...
        best.task, best.band, best.arrow, best.dz, best.p, best.q, best.N);
end

function [dz, p, q, dirArrow, Ncommon] = compute_dz_p_q(fname, sheetA, sheetB, bandCols)
% Returns 1x5 dz, p, q, arrows, and Ncommon (minimum paired N across bands)

    A = readmatrix(fname, 'Sheet', sheetA);
    B = readmatrix(fname, 'Sheet', sheetB);

    dz = nan(1, numel(bandCols));
    p  = nan(1, numel(bandCols));
    dirArrow = repmat("?", 1, numel(bandCols));
    Ns = nan(1, numel(bandCols));

    for b = 1:numel(bandCols)
        [uA, mA] = participant_means(A, bandCols(b));
        [uB, mB] = participant_means(B, bandCols(b));

        [u, ia, ib] = intersect(uA, uB);
        Ns(b) = numel(u);

        if Ns(b) < 2
            continue;
        end

        d = mB(ib) - mA(ia); % B - A

        md = mean(d, 'omitnan');
        sd = std(d, 0, 'omitnan');
        if sd ~= 0 && ~isnan(sd)
            dz(b) = md / sd;
        end

        [~, p(b)] = ttest(d, 0, 'Tail','both');
        dirArrow(b) = ternary(md >= 0, "↑", "↓");
    end

    Ncommon = min(Ns, [], 'omitnan');
    q = fdr_bh(p); % FDR across 5 bands within this task-region
end

function [uIDs, means] = participant_means(M, valueCol)
    ids = M(:,1);
    uIDs = unique(ids);
    means = nan(size(uIDs));
    for k = 1:numel(uIDs)
        means(k) = mean(M(ids==uIDs(k), valueCol), 'omitnan');
    end
end

function q = fdr_bh(p)
% Benjamini–Hochberg q-values for a p vector, returned in original order.
    p = p(:);
    q = nan(size(p));
    valid = ~isnan(p);
    pv = p(valid);
    m = numel(pv);
    if m == 0, return; end

    [ps, idx] = sort(pv, 'ascend');
    qtemp = ps .* m ./ (1:m)';
    qtemp = flipud(cummin(flipud(qtemp)));
    qtemp(qtemp > 1) = 1;

    qv = nan(size(pv));
    qv(idx) = qtemp;
    q(valid) = qv;
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
