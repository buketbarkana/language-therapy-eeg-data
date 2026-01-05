%% Stand-alone script to generate Table 4 (Task × Region therapy effects with p-values)
% Table 4. Task × Region therapy effects with p-values (TD vs DLD)
%
% For each Task–Region condition:
%   1) Compute paired differences per participant for each band: d_i = (B - A)
%   2) Compute dz = mean(d) / std(d)
%   3) Compute paired t-test p-value (two-sided) across participants
%   4) Compute q-values via FDR (Benjamini–Hochberg) across the 5 bands
%      WITHIN that Task–Region condition
%   5) Select "Top band" = band with largest |dz| (ties broken by smaller p)
%   6) Report: band (↑/↓), dz, p (q) for TD and DLD
%
% Assumptions (edit in USER SETTINGS):
%   - 12 stacked Excel files exist in dataFolder (fileList)
%   - Each file contains two sheets per group (TD and DLD) for Therapy A/B
%   - Columns: Col1=participant ID; Col2=word_id (ignored)
%   - Band power columns: delta=9, theta=10, alpha=11, beta=12, broadband=13
%
% Output:
%   - MATLAB table T4
%   - Optional Excel export

clear; clc;

%% ===== USER SETTINGS =====
dataFolder = pwd;  % change if needed

fileList = { ...
    'FWP_Frnt_stacked.xlsx','FWP_Occ_stacked.xlsx','FWP_Par_stacked.xlsx','FWP_Tmp_stacked.xlsx', ...
    'UWE_Frnt_stacked.xlsx','UWE_Occ_stacked.xlsx','UWE_Par_stacked.xlsx','UWE_Tmp_stacked.xlsx', ...
    'UWR_Frnt_stacked.xlsx','UWR_Occ_stacked.xlsx','UWR_Par_stacked.xlsx','UWR_Tmp_stacked.xlsx'};

taskOrder   = {'FWP','UWE','UWR'};
regionOrder = {'Frnt','Occ','Par','Tmp'};

bandCols  = [9 10 11 12 13];
bandNames = {'delta','theta','alpha','beta','broadband'};

% Sheet names (EDIT to match your Excel)
TD_sheetA  = "TD_TherapyA";
TD_sheetB  = "TD_TherapyB";
DLD_sheetA = "DLD_TherapyA";
DLD_sheetB = "DLD_TherapyB";

% Reporting format
digits_dz = 3;
digits_p  = 3;

% Export
exportExcel = true;
outXlsx = "Table4_TaskRegion_TopBand_dz_p_q.xlsx";
%% =========================

%% Build Task–Region keys in desired row order (12 rows)
rowLabels = cell(numel(taskOrder)*numel(regionOrder),1);
rowKeys   = cell(size(rowLabels)); % "FWP_Frnt"
k = 0;
for iT = 1:numel(taskOrder)
    for iR = 1:numel(regionOrder)
        k = k + 1;
        rowLabels{k} = sprintf('%s–%s', taskOrder{iT}, regionOrder{iR}); % en-dash
        rowKeys{k}   = sprintf('%s_%s', taskOrder{iT}, regionOrder{iR});
    end
end

%% Create file lookup
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

%% Compute Table 4 rows
T4 = table;
T4.TaskRegion = rowLabels;

TD_col  = strings(numel(rowKeys),1);
DLD_col = strings(numel(rowKeys),1);

for iRow = 1:numel(rowKeys)
    key = rowKeys{iRow};

    if ~isKey(lut, key) || ~isfile(lut(key))
        warning("Missing file for %s", key);
        TD_col(iRow)  = "NA";
        DLD_col(iRow) = "NA";
        continue;
    end

    fname = lut(key);

    % Compute per-band dz/p/q for each group
    [dz_TD,  p_TD,  q_TD,  dir_TD,  nTD ] = compute_dz_p_q(fname, TD_sheetA,  TD_sheetB,  bandCols);
    [dz_DLD, p_DLD, q_DLD, dir_DLD, nDLD] = compute_dz_p_q(fname, DLD_sheetA, DLD_sheetB, bandCols);

    % Pick top band (largest |dz|; tie -> smaller p)
    top_TD  = pick_top_band(dz_TD,  p_TD);
    top_DLD = pick_top_band(dz_DLD, p_DLD);

    % Format strings for table
    TD_col(iRow)  = format_cell(bandNames{top_TD}, dir_TD(top_TD), dz_TD(top_TD), p_TD(top_TD), q_TD(top_TD), digits_dz, digits_p, nTD);
    DLD_col(iRow) = format_cell(bandNames{top_DLD}, dir_DLD(top_DLD), dz_DLD(top_DLD), p_DLD(top_DLD), q_DLD(top_DLD), digits_dz, digits_p, nDLD);
end

T4.TD = TD_col;
T4.DLD = DLD_col;

disp("Table 4. Task × Region therapy effects with p-values (TD vs DLD)");
disp(T4);

if exportExcel
    writetable(T4, outXlsx, 'FileType','spreadsheet');
    fprintf("Saved: %s\n", outXlsx);
end

%% ===================== FUNCTIONS =====================

function [dz, p, q, dirArrow, N] = compute_dz_p_q(fname, sheetA, sheetB, bandCols)
% For a single Task-Region file and one group:
% returns 1x5 dz, p, q, arrows, and paired N

    A = readmatrix(fname, 'Sheet', sheetA);
    B = readmatrix(fname, 'Sheet', sheetB);

    dz = nan(1, numel(bandCols));
    p  = nan(1, numel(bandCols));
    dirArrow = repmat("?", 1, numel(bandCols));

    % Compute per-band paired vectors (per-participant means)
    for b = 1:numel(bandCols)
        [uA, mA] = participant_means(A, bandCols(b));
        [uB, mB] = participant_means(B, bandCols(b));

        [u, ia, ib] = intersect(uA, uB);
        N = numel(u);

        if N < 2
            dz(b) = NaN;
            p(b)  = NaN;
            dirArrow(b) = "?";
            continue;
        end

        d = mB(ib) - mA(ia); % B - A

        md = mean(d, 'omitnan');
        sd = std(d, 0, 'omitnan');

        if sd == 0 || isnan(sd)
            dz(b) = NaN;
        else
            dz(b) = md / sd;
        end

        % Two-sided paired t-test
        [~, p(b)] = ttest(d, 0, 'Tail','both');

        dirArrow(b) = ternary(md >= 0, "↑", "↓");
    end

    % FDR correction across the 5 bands within this task-region condition
    q = fdr_bh(p);
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

function idx = pick_top_band(dz, p)
% Choose band with largest |dz|; tie-break with smaller p
    absdz = abs(dz);
    if all(isnan(absdz))
        idx = 1; % default
        return;
    end
    maxVal = max(absdz, [], 'omitnan');
    cand = find(absdz == maxVal);

    if numel(cand) == 1
        idx = cand;
        return;
    end

    % tie-break: smallest p among candidates
    [~, j] = min(p(cand));
    idx = cand(j);
end

function s = format_cell(band, arrow, dz, p, q, digits_dz, digits_p, N)
% Produce: "alpha (↑), dz=0.706, p=0.086 (q=0.199)"
    if isnan(dz) || isnan(p) || isnan(q)
        s = sprintf("%s (%s), dz=NA, p=NA (q=NA) [N=%d]", band, arrow, N);
        return;
    end
    s = sprintf("%s (%s), dz=%.*f, p=%.*f (q=%.*f) [N=%d]", ...
        band, arrow, digits_dz, dz, digits_p, p, digits_p, q, N);
end

function q = fdr_bh(p)
% Benjamini-Hochberg FDR-adjusted q-values for a vector p
% Returns q in the original order.
    p = p(:);
    q = nan(size(p));

    valid = ~isnan(p);
    pv = p(valid);
    m = numel(pv);
    if m == 0
        return;
    end

    [ps, idx] = sort(pv, 'ascend');
    qtemp = ps .* m ./ (1:m)';                 % BH
    qtemp = flipud(cummin(flipud(qtemp)));     % enforce monotonicity
    qtemp(qtemp > 1) = 1;

    qv = nan(size(pv));
    qv(idx) = qtemp;

    q(valid) = qv;
end

function out = ternary(cond, a, b)
% simple ternary helper
    if cond
        out = a;
    else
        out = b;
    end
end
