%% ========== Plot each subject and mean across subjects, all days ==========
clear; clc

% in Exp3
% subj25(200.962) 1
% subj26(200.982) 2
% subj24(200.961) 3
% subj6(200.925)  4
% subj14(200.936) 5
% subj12(200.932) 6
% subj18(200.943) 7

%% ===== Parameters =====
numofsub   = 7;    % number of subjects
numofday   = 10;   % number of days
numoftr    = 60;   % trials per block
lim_RT     = 1200; % ms
lim_initDir = 100; % deg

binX = 12;% 36=10deg, 18=20deg, 12=30deg, 8=45deg, 6=60deg
binY = 6;
mw = 36; % 72=5deg, 36=10deg, 24=15deg, for mode 21

% --- colors ---
colSub = [0.75 0.75 0.75];
colDEA = [1 0 1];
colDEB = [0 0 1];
col = [192/255,192/255,192/255];


%% ===== Day 1-4, 6-9 =====

load Exp3_Training_slim.mat

% --- Plot ranges ---
rng1 = {1:10, 12:21, 23:32, 34:43};       % DEA
rng2 = {81:90, 92:101, 103:112, 114:123}; % DEB

% --- Data containers ---
dAll.DEA.RT = [];
dAll.DEB.RT = [];
dAll.DEA.initDir = [];
dAll.DEB.initDir = [];

% --- Extract metrics for each subject ---
for sub = 1:numofsub
    [DEA_RT,  DEB_RT]  = extractMetric(d, sub, numofday, numoftr, 'RT', false);
    [DEA_Dir, DEB_Dir] = extractMetric(d, sub, numofday, numoftr, 'initDir', true);

    dAll.DEA.RT      = [dAll.DEA.RT; DEA_RT];
    dAll.DEB.RT      = [dAll.DEB.RT; DEB_RT];
    dAll.DEA.initDir = [dAll.DEA.initDir; DEA_Dir];
    dAll.DEB.initDir = [dAll.DEB.initDir; DEB_Dir];
end

% === Plot ===
% --- RT ---
plotMetric(rng1, rng2, dAll.DEA.RT, dAll.DEB.RT, colSub, colDEA, colDEB, ...
           lim_RT, 'Reaction Time (ms)', 1);

% --- InitDir ---
plotMetric(rng1, rng2, dAll.DEA.initDir, dAll.DEB.initDir, colSub, colDEA, colDEB, ...
           lim_initDir, 'Initial Direction Error (deg)', 2);

% === Statistics ===
fprintf('----- RT, last 5 blocks on Day4 -----\n');
tempA = [nanmean(dAll.DEA.RT(1,36:40)),nanmean(dAll.DEA.RT(2,36:40)),nanmean(dAll.DEA.RT(3,36:40)),nanmean(dAll.DEA.RT(4,36:40)),nanmean(dAll.DEA.RT(5,36:40)),nanmean(dAll.DEA.RT(6,36:40)),nanmean(dAll.DEA.RT(7,36:40))];
tempB = [nanmean(dAll.DEB.RT(1,36:40)),nanmean(dAll.DEB.RT(2,36:40)),nanmean(dAll.DEB.RT(3,36:40)),nanmean(dAll.DEB.RT(4,36:40)),nanmean(dAll.DEB.RT(5,36:40)),nanmean(dAll.DEB.RT(6,36:40)),nanmean(dAll.DEB.RT(7,36:40))];
ttest(tempA,tempB,2,'paired')

fprintf('----- InitDir, last 5 blocks on Day4 -----\n');
fprintf('-----init Dir Error, last 5 block on Day4-----\n')
tempA = [nanmean(dAll.DEA.initDir(1,36:40)),nanmean(dAll.DEA.initDir(2,36:40)),nanmean(dAll.DEA.initDir(3,36:40)),nanmean(dAll.DEA.initDir(4,36:40)),nanmean(dAll.DEA.initDir(5,36:40)),nanmean(dAll.DEA.initDir(6,36:40)),nanmean(dAll.DEA.initDir(7,36:40))];
tempB = [nanmean(dAll.DEB.initDir(1,36:40)),nanmean(dAll.DEB.initDir(2,36:40)),nanmean(dAll.DEB.initDir(3,36:40)),nanmean(dAll.DEB.initDir(4,36:40)),nanmean(dAll.DEB.initDir(5,36:40)),nanmean(dAll.DEB.initDir(6,36:40)),nanmean(dAll.DEB.initDir(7,36:40))];
ttest(tempA,tempB,2,'paired')


%% ===== Day 5-15 =====

% Mapping: day | fileName | startsub | rng1 | rng2
dayMap = {
  5, 'Exp3_Day5ADay10B_slim.mat',        1,  46:48,       [];
  6, 'Exp3_Day6ADay10B_slim.mat',  1,  78:80,       125:127;
 10, 'Exp3_Day6ADay10B_slim.mat',  1,  78:80,       125:127;
 11, 'Exp3_Day11_slim.mat',         3,  1:3,         5:7;
 12, 'Exp3_Day12_slim.mat',         1,  [1:3,7:9],   [4:6,10:12];
 13, 'Exp3_Day13_slim.mat',         1,  1:2:11,      2:2:12;
 14, 'Exp3_Day14_slim.mat',         1,  [1:3,7:9],   [4:6,10:12];
 15, 'Exp3_Day15_slim.mat',         1,  1:2:11,      2:2:12
};

for day = 5:15

    % --- Get setting for the day --- 
    idx = find([dayMap{:,1}] == day, 1);
    if isempty(idx)
        warning('Cannot find the setting for Day %d', day);
        continue;
    end

    % --- Reset ---
    dPro.DEA.RT = []; dPro.DEA.initDir = [];
    dPro.DEB.RT = []; dPro.DEB.initDir = [];

    % --- Load --- 
    fileName = dayMap{idx,2};
    startsub = dayMap{idx,3};
    rng1     = dayMap{idx,4};
    rng2     = dayMap{idx,5};
    load(fileName); % Read dco

    if day < 12
        % Day 5〜11
        numofblk = 3;
        mapList = [1, 2]; % 1=DEA, 2=DEB
        assignPattern = {1:3, 1:3}; % 3 blocks for DEA, DEB
    else
        % Day 12〜
        numofblk = 12;
        mapList = 1; % Read Column 1, Both DEA and DEB is in Col1
        if ismember(day, [12,14])
            assignPattern = {[1:3, 7:9], [4:6, 10:12]};
        elseif ismember(day, [13, 15])
            assignPattern = {1:2:11, 2:2:12};
        end
    end

    % === Subject loop ===
    for sub = startsub:numofsub
        sRT.DEA = []; sRT.DEB = [];
        sDir.DEA = []; sDir.DEB = [];

        if day < 12
            % --- Column 1:DEA, Column 2:DEB ---
            for m = 1:2
                map = mapList(m);
                bRT = nan(1, numofblk);
                bDir = nan(1, numofblk);

                if map <= size(dco,2) && ~isempty(dco{sub,map})
                    for i = 1:numofblk
                        idxBlk = (1:60) + 60*(i-1);
                        bRT(i)  = nanmean(dco{sub,map}.RT(idxBlk));
                        bDir(i) = nanmean(abs(dco{sub,map}.initDir(idxBlk)));
                    end
                end

                if m == 1
                    sRT.DEA  = [sRT.DEA  bRT(assignPattern{1})];
                    sDir.DEA = [sDir.DEA bDir(assignPattern{1})];
                else
                    sRT.DEB  = [sRT.DEB  bRT(assignPattern{2})];
                    sDir.DEB = [sDir.DEB bDir(assignPattern{2})];
                end
            end

        else
            % --- Column 1: DEA & DEB ---
            map = mapList(1);
            bRT = nan(1, numofblk);
            bDir = nan(1, numofblk);

            if map <= size(dco,2) && ~isempty(dco{sub,map})
                for i = 1:numofblk
                    idxBlk = (1:60) + 60*(i-1);
                    bRT(i)  = nanmean(dco{sub,map}.RT(idxBlk));
                    bDir(i) = nanmean(abs(dco{sub,map}.initDir(idxBlk)));
                end
            end

            % Split into patterns
            sRT.DEA  = [sRT.DEA  bRT(assignPattern{1})];
            sDir.DEA = [sDir.DEA bDir(assignPattern{1})];
            sRT.DEB  = [sRT.DEB  bRT(assignPattern{2})];
            sDir.DEB = [sDir.DEB bDir(assignPattern{2})];
        end

        % --- All participants ---
        dPro.DEA.RT       = [dPro.DEA.RT; sRT.DEA];
        dPro.DEB.RT       = [dPro.DEB.RT; sRT.DEB];
        dPro.DEA.initDir  = [dPro.DEA.initDir; sDir.DEA];
        dPro.DEB.initDir  = [dPro.DEB.initDir; sDir.DEB];
    end

    % === Plot RT ===
    plotMetric_AB(day, dPro, 'RT', col, lim_RT, rng1, rng2, 'Reaction Time (ms)', 1);

    % === Plot initDir ===
    plotMetric_AB(day, dPro, 'initDir', col, lim_initDir, rng1, rng2, 'Absolute Initial Direction Error (deg)', 2);

end

%% ===== Day 11, 16, 17 =====

% Day 11, switch by trials, short 1(testing 1)
% Day 16, switch by trials, long
% Day 17, switch by trials, short 2(testing 2)

load Exp3_Switch.mat
SW = dAll_Exp3_SW;

% --- Steady State Map A, B ---
load Exp3_Day5ADay10B_slim.mat % 1: DEA on Day5, 2: DEB on Day 10 
SSA = arrayfun(@(s) nanmean(abs(dco{s,1}.initDir)), 1:7);
SSB = arrayfun(@(s) nanmean(abs(dco{s,2}.initDir)), 1:7);
SSA_B = arrayfun(@(s) nanmean(abs(dco{s,1}.initDir_DEB)), 1:7);
SSB_A = arrayfun(@(s) nanmean(abs(dco{s,2}.initDir_DEA)), 1:7);

% --- Policy Recovery --- 
load Exp3_PolicyRecovery.mat
DEA_DEA = mean(mean(Result{1,1},2)); DEA_DEB = mean(mean(Result{1,2},2)); DEA_N = mean(mean(Result{1,3},2));
DEB_DEA = mean(mean(Result{2,1},2)); DEB_DEB = mean(mean(Result{2,2},2)); DEB_N = mean(mean(Result{2,3},2));

subjIDs = {
    '200962', '200982', '200925', '200961', '200936', '200932', '200943'
};

% Mapping:day | cond | startsub | numofsub | trialAfterSW | rng3 | rng5
dayMap = {
    11, 3, 3, 7, 10, 3:12, 16:25;
    17, 2, 1, 7, 10, 3:12, 16:25;
    16, 1, 1, 7, 15, 3:17, 20:34
};

% Chose data for map *********
SSchoice = 1; % 0: Day 5/10(co-task under BA/DE), 1: Day 16(> 6 trials post switch)

for day = [11, 16, 17]

    % ===== Initial Direction Error and RT =====
    % --- Get setting for the day ---
    idx = find([dayMap{:,1}] == day, 1);
    if isempty(idx)
        warning('Cannot find the setting for Day %d', day);
        continue;
    end

    cond            = dayMap{idx,2};
    startsub        = dayMap{idx,3};
    numofsub        = dayMap{idx,4};
    trialAfterSW    = dayMap{idx,5};
    rng3            = dayMap{idx,6};
    rng5            = dayMap{idx,7};

    % --- All participants ---
    dBtoA.RT      = [];
    dAtoB.RT      = [];
    dBtoA.initDir = [];
    dAtoB.initDir = [];
    for sub = startsub:numofsub
        SW{cond,sub}.BtoA.initDir = abs(SW{cond,sub}.BtoA.initDir);
        SW{cond,sub}.AtoB.initDir = abs(SW{cond,sub}.AtoB.initDir);

        dBtoA.RT      = [dBtoA.RT;      SW{cond,sub}.BtoA.RT];
        dAtoB.RT      = [dAtoB.RT;      SW{cond,sub}.AtoB.RT];
        dBtoA.initDir = [dBtoA.initDir; SW{cond,sub}.BtoA.initDir];
        dAtoB.initDir = [dAtoB.initDir; SW{cond,sub}.AtoB.initDir];
    end

    % --- Plot RT ---
    plotMetric_SwitchUnified(day, ...
        dBtoA, dAtoB, [], [],[], [], ...        % data (switching, steady state)
        'RT', col, lim_RT, ...   % metrics, color, range of Y axis
        rng3, rng5, ...    % range of X axsis
        'Reaction Time (ms)', 1) % label, Figure No.

    % --- Plot Init Dir Error --- 
    plotMetric_SwitchUnified(day, ...
        dBtoA, dAtoB, SSA, SSB, SSA_B, SSB_A, ...        
        'initDir', col, lim_initDir, ...   
        rng3, rng5, ...    
        'Absolute Initial Direction Error (deg)', 2)  

    % ===== Policy =====
    % --- Initialization ---
    for sw = 1:2
        pAB_opt{sw} = [];
        eMix{sw} = [];
        for p = 1:3
            Result{sw,p} = [];
        end
        for p = 1:2
            rePol{sw,p} = [];
        end
    end

    for sub = startsub:numofsub

        fprintf('Subject %d: %s\n', sub, subjIDs{sub});

        [targAng, initDir] = load_steady_maps(sub, SW, SSchoice, trialAfterSW, day);

        [targAng, initDir] = wrap_angles(targAng, initDir, trialAfterSW);
         
        p_map = build_probability_maps(targAng, initDir, mw, binX, binY, trialAfterSW);

        pAB_opt = estimate_trial_policies(SW, cond, sub, p_map, trialAfterSW);

        % --- All participants ---
        for sw = 1:2 % 1=BtoA, 2=AtoB
            Result{sw,1}(sub,:) = pAB_opt{sw}(:,1); % policy DEA
            Result{sw,2}(sub,:) = pAB_opt{sw}(:,2); % policy DEB
            Result{sw,3}(sub,:) = 1-(pAB_opt{sw}(:,1)+pAB_opt{sw}(:,2)); % other policy 
        end

    end

    % eMix = p(persist) * e[error|persist] + (p(correct)+p(other)) * e[error|no-persist]
    % policy: Result{:,1}: DEA, Result{:,2}: DEB
    % direction: Result{1,:}: DEB->DEA, Result{2,:}: DEA->DEB
    for sw = 1:2
        if sw == 1 % DEB to DEA
            rePol{sw,1} = Result{sw,1} + Result{sw,3}; % correct + other
            rePol{sw,2} = Result{sw,2}; % persist
            eMix{sw} = SSB_A' .* rePol{sw,2} + SSA' .* rePol{sw,1}; % SSB_A: DEB recunstructed under DEA
        else % DEA to DEB
            rePol{sw,1} = Result{sw,1}; % persist
            rePol{sw,2} = Result{sw,2} + Result{sw,3}; % correct + other
            eMix{sw} = SSA_B' .* rePol{sw,1} + SSB' .* rePol{sw,2}; % SSA_B: DEA reconstruct under DEB
        end
    end

    % ===== Group results: policy =====
    figure(2); hold on;
    for sw = 1:2 % 1=BtoA, 2=AtoB
        rngUse = ternary(sw==1, rng3, rng5);
        plotPolicyWithError(day, rngUse, Result, sw, startsub, trialAfterSW, rng3, rng5, DEA_DEA, DEA_DEB, DEA_N, DEB_DEA, DEB_DEB, DEB_N)
        ploteMix(day, sw, eMix, rng3, rng5, startsub, numofsub);
    end

end



%% ===== function - Training =====

% === Get data ===
function [DEA, DEB] = extractMetric(data, subj, numofday, numoftr, metric, absFlag)

    DEA = [];
    DEB = [];

    for day = 1:numofday
        % --- Check block number ---
        if day == 1
            numofblk = 13;
        elseif day == 5 || day == 10
            numofblk = 2;
        else
            numofblk = 11;
        end

        blockMean = nan(1, numofblk);
        for i = 1:numofblk
            vals = data{subj,day}.(metric)(1+numoftr*(i-1):numoftr*i);
            if absFlag
                vals = abs(vals);
            end
            blockMean(i) = nanmean(vals);
        end

        % Divided into DEA / DEB on each day 
        if day == 1
            DEA = [DEA blockMean(4:numofblk)];
        elseif day > 1 && day < 6
            DEA = [DEA blockMean(2:numofblk)];
        elseif day >= 6
            DEB = [DEB blockMean(2:numofblk)];
        end
    end

end
 

% === Plot Training DEA(Day1-4) / DEB(6-9) ===
function plotMetric(rng1, rng2, DEA, DEB, colorSub, colorDEA, colorDEB, limY, yLabelText, figNum)

    numofsub = size(DEA,1);

    % --- Settings for Figure ---
    fhandle = figure(figNum); clf; hold on;
    set(fhandle, 'Units', 'inch', ...
                 'Position', [100, 100, 11.69, 8.27], ...
                 'Color', 'w');

    subplot(4,10,1:8); hold on;
    for d = 1:4
        for sub = 1:numofsub
            plot(rng1{d}, DEA(sub,1+10*(d-1):10*d), '-', 'color', colorSub, 'LineWidth', 1.2);
            plot(rng2{d}, DEB(sub,1+10*(d-1):10*d), '-', 'color', colorSub, 'LineWidth', 1.2);
        end
        plot(rng1{d}, nanmean(DEA(:,1+10*(d-1):10*d)), 'o-', 'Color', colorDEA, 'MarkerFaceColor', colorDEA, 'MarkerSize', 3, 'LineWidth', 2);
        plot(rng2{d}, nanmean(DEB(:,1+10*(d-1):10*d)), 'o-', 'Color', colorDEB, 'MarkerFaceColor', colorDEB, 'MarkerSize', 3, 'LineWidth', 2);
    end

    % --- Set vertical bars ---
    vlines = [11 22 33 44 91 102 113 124];

    % --- Set axes ---
    if contains(lower(yLabelText), 'reaction') % RT
        for v = vlines
            plot([v v],[0 limY],'k-');
        end
        set_axes(1, limY, 'reaction', yLabelText)
    elseif contains(lower(yLabelText), 'initial') % InitDir
        for v = vlines
            plot([v v],[0 limY],'k-');
        end
        set_axes(1, limY, 'initial', yLabelText)
    end

end


% === Plot DEA, DEB, Practice Switching, Day 5,6,10-15 ===
function plotMetric_AB(day, dPro, metrics, col, limY, rng1, rng2, yLabelText, figNum)

    DEA_data = dPro.DEA.(metrics);
    DEB_data = dPro.DEB.(metrics);

    figure(figNum); hold on;

    switch day
        case {5, 6} % DEA only
            subplot(4,10,1:8); hold on;
            plotSubjects(DEA_data, rng1, col, '-');
            plotMean(DEA_data, rng1, 'm', 'o-', 2);

        case 10 % DEB only
            subplot(4,10,1:8); hold on;
            plotSubjects(DEB_data, rng2, col, '-');
            plotMean(DEB_data, rng2, 'b', 'o-', 2);

        case 11 % DEA & DEB
            subplot(4,10,11:12); hold on;
            plotSubjects(DEA_data, rng1, col, '-');
            plotSubjects(DEB_data, rng2, col, '-');
            plotMean(DEA_data, rng1, 'm', 'o-', 2);
            plotMean(DEB_data, rng2, 'b', 'o-', 2);
            set_axes(day, limY, metrics, yLabelText)

        case {12, 14} % Split pattern
            if day == 12
                subplot(4,10,13:14); hold on;
            else
                subplot(4,10,17:18); hold on;
            end
            plotSplitPattern(DEA_data, DEB_data, rng1, rng2, col, '-');
            plotSplitMean(DEA_data, DEB_data, rng1, rng2, 'o-');
            set_axes(day, limY, metrics, yLabelText)

        case {13, 15} % Marker only
            if day == 13
                subplot(4,10,15:16); hold on;
            else
                subplot(4,10,19:20); hold on;
            end
            plotSubjects(DEA_data, rng1, col, 'o');
            plotSubjects(DEB_data, rng2, col, 'o');
            plotMean(DEA_data, rng1, 'm', 'o', 1.5);
            plotMean(DEB_data, rng2, 'b', 'o', 1.5);
            set_axes(day, limY, metrics, yLabelText)
    end
end


% === Plot individuals ===
function plotSubjects(data, xVals, col, style)
    for sub = 1:size(data,1)
        plot(xVals, data(sub,:), style, 'Color', col, 'MarkerFaceColor', col, ...
        'MarkerSize', 3, 'LineWidth', 1.2);
    end
end


% === Plot group mean ===
function plotMean(data, xVals, col, style, lw)
    plot(xVals, nanmean(data), style, 'Color', col, 'MarkerFaceColor', col, ...
        'MarkerSize', 3, 'LineWidth', lw);
end


% === Plot individuals, switching by blocks ===
function plotSplitPattern(DEA, DEB, rng1, rng2, col, style)
    for sub = 1:size(DEA,1)
        plot(rng1(1:3), DEA(sub,1:3), style, 'Color', col, 'LineWidth', 1.2);
        plot(rng2(1:3), DEB(sub,1:3), style, 'Color', col, 'LineWidth', 1.2);
        plot(rng1(4:6), DEA(sub,4:6), style, 'Color', col, 'LineWidth', 1.2);
        plot(rng2(4:6), DEB(sub,4:6), style, 'Color', col, 'LineWidth', 1.2);
    end
end


% === Plot group mean, switching by blocks ===
function plotSplitMean(DEA, DEB, rng1, rng2, style)
    plot(rng1(1:3), nanmean(DEA(:,1:3)), style, 'Color', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', 3, 'LineWidth', 2);
    plot(rng2(1:3), nanmean(DEB(:,1:3)), style, 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 3, 'LineWidth', 2);
    plot(rng1(4:6), nanmean(DEA(:,4:6)), style, 'Color', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', 3, 'LineWidth', 2);
    plot(rng2(4:6), nanmean(DEB(:,4:6)), style, 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 3, 'LineWidth', 2);
end



%% ====== function - Switching (init dir error & RT) =====

% --- Plot init dir error & RT ---
function plotMetric_SwitchUnified(day, dBtoA, dAtoB, SSA, SSB, SSA_B, SSB_A, metrics, col, limY, rng3, rng5, yLabelText, figNum)

    dataB = dBtoA.(metrics); % B→A
    dataA = dAtoB.(metrics); % A→B

    figure(figNum); hold on;

    % --- Set subplot location ---
    switch day
        case 11
            sp = 21:23;
            trpost = 1:10; % trial post switch 
        case 16
            sp = 24:27; 
            trpost = 1:15;
        case 17
            sp = 28:30; 
            trpost = 1:10;
    end

    % --- Plot B->A, A->B ---
    subplot(4,10,sp); hold on;
    plotSubjects_switch(dataB, rng3, trpost, col, '-');
    plotMean_switch(dataB, rng3, trpost, 'm', 'o-', 2);
    plotSubjects_switch(dataA, rng5, trpost, col, '-');
    plotMean_switch(dataA, rng5, trpost, 'b', 'o-', 2);

    if strcmp(metrics, 'initDir')
        % steadystate A
        plot(rng3,ones(1,length(rng3))*nanmean(SSA),'m-','Color','m','LineWidth',1);
        plot(rng3,ones(1,length(rng3))*nanmean(SSB_A),'b:','Color','b','LineWidth',1.2);
        % steadystate B
        plot(rng5,ones(1,length(rng5))*nanmean(SSB),'b-','Color','b','LineWidth',1);
        plot(rng5,ones(1,length(rng5))*nanmean(SSA_B),'m:','Color','m','LineWidth',1.2);
    end

    set_axes(day, limY, metrics, yLabelText);

end


% === Plot individuals ===
function plotSubjects_switch(dataMartix, xVals, trange, col, style)
    for sub = 1:size(dataMartix,1)/50 % 50 tr per subject
        plot(xVals, nanmean(dataMartix(1+50*(sub-1):50*sub,trange)), style, 'Color', col, 'LineWidth', 1.2);
    end
end


% === Plot group mean ===
function plotMean_switch(dataMatrix, xVals, trange, col, style, lw)
    plot(xVals, nanmean(dataMatrix(:, trange)), style, 'Color', col, 'MarkerFaceColor', col, ...
        'MarkerSize', 3, 'LineWidth', lw);
end



%% ===== Set axes =====

function set_axes(day, limY, metrics, yLabelText)
    ylim([0 limY])
    if any( contains(string(metrics), ["reaction","RT"], 'IgnoreCase', true) ) % RT
        yticks([0 300 600 900 1200])
        if day == 1 || day == 11
            yticklabels([0 300 600 900 1200])
            ylabel(yLabelText);
        else
            yticklabels([]);
        end
    else % init dir error
        yticks([0 25 50 75 100])
        if day == 1 || day == 11
            yticklabels([0 25 50 75 100])
            ylabel(yLabelText);
        else
            yticklabels([]);
        end
    end
    if day == 1
        xticks([])
        xticklabels([])
        xlabel('Day');
    else 
        xticks(1:34)
        xticklabels([])
        xlabel('Block');
    end
end


%% ===== function - Switching (policy) =====

% === Get distribution from center-out task ===
function [targAng, initDir] = load_steady_maps(sub, SW, SSchoice, trialAfterSW, day)

    switch SSchoice
        case {0} % Day 5/10 for map
            load Exp3_Day5ADay10B_slim.mat
            for sw = 1:2
                for tr = 1:trialAfterSW
                    targAng{tr}{sw} = dco{sub,sw}.targAng;
                    initDir{tr}{sw,1} = dco{sub,sw}.initDirL; % hd = 1:2 % 1:LH, 2:RH
                    initDir{tr}{sw,2} = dco{sub,sw}.initDirR;
                end
            end
    
        case {1} % Day 16 for map
            maps = {'BtoA', 'AtoB'};
    
            for k = 1:numel(maps)
                temp1{k}    = SW{1,sub}.(maps{k}).targAng; % SW{1,:}->Day16
                temp2{k,1}  = SW{1,sub}.(maps{k}).initDirL;
                temp2{k,2}  = SW{1,sub}.(maps{k}).initDirR;
            end
      
            switch day
                case {16} % DEA only
                    for sw = 1:2
                        for tr = 1:trialAfterSW
                            if tr < 6
                                cols = 6:15;
                            else
                                cols = 6:15;
                                cols(tr-5) = []; % delete row used for validation
                            end

                            % targAng
                            targAng{tr}{sw} = reshape(temp1{sw}(:, cols), 1, []);

                            % initDir
                            for hd = 1:2  % 1:LH, 2:RH
                                initDir{tr}{sw,hd} = reshape(temp2{sw,hd}(:, cols), 1, []);
                            end
                        end
                    end
                        
                case {11, 17} % use 6-15 trials post switch for map
                    cols = 6:15;
                    for sw = 1:2
                        for tr = 1:trialAfterSW
                            % targAng
                            targAng{tr}{sw} = reshape(temp1{sw}(:, cols), 1, []);

                            % initDir
                            for hd = 1:2  % 1:LH, 2:RH
                                initDir{tr}{sw,hd} = reshape(temp2{sw,hd}(:, cols), 1, []);
                            end
                        end
                    end
            end
  
    end
end


% === Edged wraped up ===
% add points to the edges and use different bin size for x and y
% -200 to 200 deg
function [targAng, initDir] = wrap_angles(targAng, initDir, trialAfterSW)
    for sw = 1:2
        for tr = 1:trialAfterSW
            idx_low  = targAng{tr}{sw} > -180 & targAng{tr}{sw} < -160;
            idx_high = targAng{tr}{sw} > 160  & targAng{tr}{sw} < 180;
    
            % wrap-around target angles
            extraTarg{tr} = [targAng{tr}{sw}(idx_low) + 360, targAng{tr}{sw}(idx_high) - 360];
    
            % wrap-around initDir for both hands
            for hd = 1:2
                extraInit{tr}{hd} = [initDir{tr}{sw,hd}(idx_low), initDir{tr}{sw,hd}(idx_high)];
                initDir{tr}{sw,hd} = [initDir{tr}{sw,hd}, extraInit{tr}{hd}];
            end
    
            targAng{tr}{sw} = [targAng{tr}{sw}, extraTarg{tr}];
        end
    end
end


% === Build probability map ===
function p_map = build_probability_maps(targAng, initDir, mw, binX, binY, trialAfterSW)

    % --- Smooting (get moving average) ---
    for sw = 1:2
        for tr = 1:trialAfterSW
            for hd = 1:2 % hand, 1:left hand, 2:right hand
                HtargMat{tr}{sw,hd} = zeros(100,100);
            end
            for hd = 1:2
                for i = 1:size(targAng{tr}{sw},2)
                    cntT = 1;
                    while -200 + (cntT-1)*360/mw + 360/binX <= 200
                        if  targAng{tr}{sw}(i)>-200+(cntT-1)*360/mw && targAng{tr}{sw}(i)<=-200+(cntT-1)*360/mw+360/binX && ~isnan(targAng{tr}{sw}(i))
                            cntH = 1;
                            while -200 + (cntH-1)*360/mw + 360/binY <= 200
                                if initDir{tr}{sw,hd}(i)>-200+(cntH-1)*360/mw && initDir{tr}{sw,hd}(i)<=-200+(cntH-1)*360/mw+360/binY && ~isnan(initDir{tr}{sw,hd}(i))
                                    HtargMat{tr}{sw,hd}(cntH,cntT) = HtargMat{tr}{sw,hd}(cntH,cntT) + 1;
                                end
                                cntH = cntH + 1;
                            end
                        end
                        cntT = cntT + 1;
                    end
                end
            end
        end
    end
    cntT = cntT - 1;
    cntH = cntH - 1;
    
    % --- Normalization ---
    for sw = 1:2
        for tr = 1:trialAfterSW
            for hd = 1:2
                colSum = sum(HtargMat{tr}{sw,hd}, 1);
                ProHtargMat{tr}{sw,hd} = HtargMat{tr}{sw,hd} ./ (colSum + (colSum==0));  % avoid div/0
            end
        end
    end
    
    % --- Final policy maps ---
    for tr = 1:trialAfterSW
        p_map{tr}{1,1} = ProHtargMat{tr}{1,1}(1:cntH,1:cntT); % mapA, left
        p_map{tr}{1,2} = ProHtargMat{tr}{1,2}(1:cntH,1:cntT); % mapA, right
        p_map{tr}{2,1} = ProHtargMat{tr}{2,1}(1:cntH,1:cntT); % mapB, left
        p_map{tr}{2,2} = ProHtargMat{tr}{2,2}(1:cntH,1:cntT); % mapB, right
    end

end


% === Estimate maximum liklihood, optimize pA, pB ===
function pAB_opt = estimate_trial_policies(SW, cond, sub, p_map, trialAfterSW)

    for sw = 1:2 %1:mapB to mapA, 2:mapA to mapB
    
        % --- theta: all trials post switch ---
        for tr = 1:trialAfterSW

            % L:1, R:2, T:3
            % HtargMat{sw,hd}
            if sw == 1
                theta = [];
                theta(:,1) = SW{cond,sub}.BtoA.initDirL(:,tr);
                theta(:,2) = SW{cond,sub}.BtoA.initDirR(:,tr);
                theta(:,3) = SW{cond,sub}.BtoA.targAng(:,tr);
            elseif sw == 2
                theta = [];
                theta(:,1) = SW{cond,sub}.AtoB.initDirL(:,tr);
                theta(:,2) = SW{cond,sub}.AtoB.initDirR(:,tr);
                theta(:,3) = SW{cond,sub}.AtoB.targAng(:,tr);
            end

            % --- Find rows where any column has 0 or NaN ---
            % 0: one of hands doesn't move
            rowsToDelete = any(theta == 0 | isnan(theta), 2);
            theta(rowsToDelete, :) = []; % Remove those rows
    
            % --- Define binning parameters ---
            binSizeT = 10; offsetT = -185;  % Target angle bins: from -185 to 200
            binSizeH = 10; offsetH = -170;  % Hand angle bins: from -170 to 200
    
            % --- Compute bin indices ---
            I = zeros(size(theta));  % [N×3] matrix for bin indices
    
            % --- Convert angles to bin indices ---
            % count trials in each bin (target, right hand, left hand)
            I(:,3) = floor((theta(:,3) - offsetT) / binSizeT) + 1;  % Target
            I(:,1) = floor((theta(:,1) - offsetH) / binSizeH) + 1;  % Left hand
            I(:,2) = floor((theta(:,2) - offsetH) / binSizeH) + 1;  % Right hand
    
            % --- Validate indices ---
            valid = all(I > 0 & mod(I,1) == 0, 2);
            I = I(valid, :);
            theta = theta(valid, :);
    
            % --- Compute number of bins ---
            CNTT = floor((200 - offsetT) / binSizeT);  % Number of target bins
            CNTH = floor((200 - offsetH) / binSizeH);  % Number of hand bins
    
            % --- Compute liklihood for each trial ---
            p_I = zeros(size(theta,1),2);
            for n = 1:size(theta,1)
                p_I(n,1) = p_map{tr}{1,1}(I(n,1),I(n,3)) * p_map{tr}{1,2}(I(n,2),I(n,3)) - (1/CNTT)^2;
                p_I(n,2) = p_map{tr}{2,1}(I(n,1),I(n,3)) * p_map{tr}{2,2}(I(n,2),I(n,3)) - (1/CNTT)^2;
                % p_I(n,3) = (1/360) * (1/360);
            end
    
            % --- Maximum liklihood estimation：optimize pA, pB ---
            pAB = [.45; .45]; % [pA; pB]
            A = [1,1]; b = 1; lb = [1e-6,1e-6]; ub = [1,1];
            LL = @(params) -sum(log(p_I*params + (1/CNTT)^2));
            LLv = @(params) -log(p_I*params + (1/CNTT)^2);
            pAB_opt{sw}(tr,:) = fmincon(LL,pAB,A,b,[],[],lb,ub);
    
        end
    end
end


% === Plot probability of policy ===
function plotPolicyWithError(day, rngVals, Result, sw, startIdx, trialAfterSW, rng3, rng5, DEA_DEA, DEA_DEB, DEA_N, DEB_DEA, DEB_DEB, DEB_N)

    % --- Set subplot ---
    switch day
        case 11
            sp = 31:33;
        case 16
            sp = 34:37;
        case 17
            sp = 38:40;
    end

    subplot(4,10,sp); hold on;
    if day == 11
        ylabel('Probability');
    end
    xlabel('Trial Post Switch');
    xticks(1:34)
    xticklabels([])
    yticks(0:0.2:1);
    yticklabels(0:0.2:1);

    % --- Plot policy -----
    colors = {'g','m','b'}; % other, de novo, base
    order  = [3,1,2];       % Result{sw,3}, Result{sw,1}, Result{sw,2}
    for i = 1:length(order)
        shadedErrorBar( ...
            rngVals, ...
            nanmean(Result{sw,order(i)}(startIdx:end,1:trialAfterSW)), ...
            seNaN(Result{sw,order(i)}(startIdx:end,1:trialAfterSW)), ...
            {sprintf('%s-o',colors{i}), ...
            'markerfacecolor',colors{i}, ...
            'markeredgecolor',colors{i}, ...
            'markersize',3, 'LineWidth',2} ...
        );
    end

    % --- Plot policy recovery ---
    if sw == 1 % DEBtoDEA
        rngUse = rng3;
        plotPolicyRecovery(rngUse, DEA_DEA, DEA_DEB, DEA_N);
    else
        rngUse = rng5;
        plotPolicyRecovery(rngUse, DEB_DEA, DEB_DEB, DEB_N);
    end
 
end


% === Plot policy recovery ===
function plotPolicyRecovery(rngVals, valDE, valBA, valOther)
    plot(rngVals, ones(size(rngVals))*valDE,    'm--','LineWidth',1); % de novo
    plot(rngVals, ones(size(rngVals))*valBA,    'b--','LineWidth',1); % base
    plot(rngVals, ones(size(rngVals))*valOther, 'g--','LineWidth',1); % other
end


% === Plot eMix ===
function ploteMix(day, sw, eMix, rng3, rng5, startsub, numofsub)

    % --- Set subplot location ---
    switch day
        case 11
            sp = 21:23;
        case 16
            sp = 24:27; 
        case 17
            sp = 28:30; 
    end
    subplot(4,10,sp); hold on;
    
    % --- Plot eMix ---
    if sw == 1
        plot(rng3, nanmean(eMix{sw}(startsub:numofsub,:)), 'ro-','LineWidth',1.2);
    else
        plot(rng5, nanmean(eMix{sw}(startsub:numofsub,:)), 'ro-','LineWidth',1.2);
    end
    

end

function out = ternary(cond, valTrue, valFalse)
    if cond, out = valTrue; else, out = valFalse; end
end

return;


