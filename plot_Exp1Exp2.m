%% Choice of experiment you want to plot

clear; clc
expID = 2; % 1:Exp 1, 2:Exp 2

%% Exp 1 (training)

if expID == 1
    load Exp1_Training.mat

    numofsub = 25;
    numofday = 4;
    numoftr = 60;

    dAll.BA.initDir = [];
    dAll.DE.initDir = [];

    xFigPos = 1200; yFigPos = 200; figWidth = 900/2; figHeight = 565/2;
    figure(10); clf; set(gcf,'Color', 'w', 'Units', 'points', 'Position', [xFigPos yFigPos figWidth figHeight]);
    hold on;
    col = [192/255,192/255,192/255];
    lim = 90;

    for subj = 1:numofsub

        bDir = zeros(1,10);
        sDir.BA = []; sDir.DE = [];

        for day = 1:numofday

            if day == 1
                numofblk = 13;
            else
                numofblk = 11;
            end

            %%% init dir
            for i = 1:numofblk
                bDir(1,i) = nanmean(abs(d{subj,day}.initDir(1+60*(i-1):60*i))); % absolute mean
            end
            if day == 1
                sDir.BA = [sDir.BA bDir(1,1:3)];
                sDir.DE = [sDir.DE bDir(1,4:numofblk)];
            else
                sDir.BA = [sDir.BA bDir(1,1)];
                sDir.DE = [sDir.DE bDir(1,2:numofblk)];
            end


        end

        % all subjects data
        dAll.BA.initDir = [dAll.BA.initDir; sDir.BA];
        dAll.DE.initDir = [dAll.DE.initDir; sDir.DE];

    end

    %%% initDir
    for d = 1:4
        for sub = 1:numofsub
            plot(1+10*(d-1):10*d, dAll.DE.initDir(sub,1+10*(d-1):10*d),'-','color',col,'LineWidth',1);
        end
        shadedErrorBar(1+10*(d-1):10*d,nanmean(dAll.DE.initDir(:,1+10*(d-1):10*d))',seNaN(dAll.DE.initDir(:,1+10*(d-1):10*d))','m.-',3);
    end
    shadedErrorBar((1:3+(numofday-1))+numofday*10,nanmean(dAll.BA.initDir(:,1:3+(numofday-1)))',seNaN(dAll.BA.initDir(:,1:3+(numofday-1)))','b.-',3); % base
    plot(repmat((10.5:10:10*(numofday+1))',1,2),[-lim lim],'k')
    plot(repmat((5:10:10*(numofday+1)-5)',1,2),[-lim lim],'k:')

    axis([0 10*numofday+10 0 lim])
    xticks([10.5 20.5 30.5 40.5])
    xticklabels([])
    yticks([0 30 60 90])
    xlabel('Block')
    ylabel('Initial Direction Error (deg)')

end


%% Exp 1 (Switching test)
% mapB(baseline) and mapA(de novo)

[cfg, totSW, SS, reps, numSubjects, dataCO, PR] = load_experiment(expID);

% === Parameters ===
trialAfterSW = 10;
directions = {'BtoA','AtoB'};  % Switch directions, BasetoDenovo, DenovotoBase
metrics    = {'initDir','RT'}; % Metrics to plot

binX = 12; % 36=10deg, 18=20deg, 12=30deg, 8=45deg, 6=60deg
binY = 6;
mw   = 36; % 72=5deg, 36=10deg, 24=15deg, for mode 21


% === Plot Settings ===
col        = [192 192 192] / 255; % subject line color
colors     = struct('BtoA','m','AtoB','b'); % Line colors
steadyIdx  = struct('BtoA',2,'AtoB',1);     % Index for steady state in D5ASS
xRanges    = struct('initDir',[0 30],'RT',[0 13]); % X-axis ranges
yRanges    = struct('initDir',[0 120],'RT',[200 800]); % Y-axis ranges

xFigPos = -30.0; yFigPos = 10.0; xFigSize = 8.27; yFigSize = 11.69;

xtickVals = struct( ...
    'initDir', [zeros(1,3),1,zeros(1,3),5,zeros(1,4),10,zeros(1,7),1,zeros(1,3),5,zeros(1,4),10], ...
    'RT',      [zeros(1,3),1,zeros(1,3),5,zeros(1,4),10] );

subplotMap = struct( ...
    'initDir', [18:23; 18:23], ... % dummy, will map below
    'RT',      [18:20; 22:24] );

trialRanges.initDir = {3:12, 20:29};
trialRanges.RT      = {3:12, 3:12};

steadyX = struct( ...
    'initDir', [2.5 13; 19.5 30], ...
    'RT',      [2.5 13; 2.5 13] );

subplotPos.initDir = { 2:7, 10:15 };
subplotPos.RT_BtoA = { 2:4, 10:12 };
subplotPos.RT_AtoB = { 6:8, 14:16 };

% === Create figures ===
for i = 1:3
    figure(i); clf; set(gcf,'units','inch','position',[xFigPos,yFigPos,xFigSize,yFigSize],'Color','w'); % InitDir Figure
end


%% ===== Combine all subjects' data =====
% --- Initialize ---
for d = directions
    for m = metrics
        dCombined.(m{1}).(d{1}) = [];
        mPerSubj.(m{1}).(d{1})  = [];
    end
end

for subj = 1:numSubjects
    % --- Ensure initDir is absolute ---
    for d = directions
        totSW{subj}.(d{1}).initDir = abs(totSW{subj}.(d{1}).initDir);
    end
    
    % --- Concatenate all trials across subjects ---
    for d = directions
        for m = metrics
            dCombined.(m{1}).(d{1}) = [dCombined.(m{1}).(d{1}); totSW{subj}.(d{1}).(m{1})];
        end
    end
    
    % --- Mean per subject (all trials) ---
    for d = directions
        for m = metrics
            mPerSubj.(m{1}).(d{1}) = [mPerSubj.(m{1}).(d{1}); nanmean(totSW{subj}.(d{1}).(m{1}))];
        end
    end
end


%% ===== Probability of Policy ===== 

% --- Initialization ---
for sw = 1:2
    pAB_opt{sw} = [];
    for p = 1:3
        Result{sw,p} = [];
    end
end

for sub = 1:numSubjects

    fprintf('Subject %d\n', sub);

    [targAng, initDir] = load_steady_maps(sub, dataCO);

    [targAng, initDir] = wrap_angles(targAng, initDir);

    p_map = build_probability_maps(targAng, initDir, mw, binX, binY);
    p_map_all{sub} = p_map;

    pAB_opt = estimate_trial_policies(totSW, sub, p_map, trialAfterSW);

    % --- All participants ---
    for sw = 1:2 % 1=BtoA, 2=AtoB
        Result{sw,1}(sub,:) = pAB_opt{sw}(:,1); % policy De Novo
        Result{sw,2}(sub,:) = pAB_opt{sw}(:,2); % policy Baseline
        Result{sw,3}(sub,:) = 1-(pAB_opt{sw}(:,1)+pAB_opt{sw}(:,2)); % other policy
    end

end

% eMix = p(persist) * e[error|persist] + (p(correct)+p(other)) * e[error|no-persist]
% policy: Result{:,1}: Denovo, Result{:,2}: Baseline
% direction: Result{1,:}: Base->Denovo, Result{2,:}: Denovo->Base
for d = 1:length(directions)
    switch d 
        case 1 % Base->Denovo, BtoA
            rePol{d,1} = Result{d,1} + Result{d,3};
            rePol{d,2} = Result{d,2};
            eMix{d} = rePol{d,2} .* SS.initDir_DE(:,1) + rePol{d,1} .* SS.initDir(:,2); % SS.initDir_DE(:,1): BA reconstructed under DE
        case 2 % Denovo->Base, AtoB
            rePol{d,1} = Result{d,1};
            rePol{d,2} = Result{d,2} + Result{d,3};
            eMix{d} = rePol{d,1} .* SS.initDir_BA(:,2) + rePol{d,2} .* SS.initDir(:,1); % SS.initDir_BA(:,2): DE reconstructed under BA
    end
end



%% ===== Plot Group =====
for m = 1:length(metrics) % initDir / RT
    for d = 1:length(directions) % BaseToDenovo / DenovoToBase
        
        % --- Set figure handle, location, x axsis etc for plot ---
        if strcmp(metrics{m}, 'initDir')
            figure(1);
            subplot(4,8, subplotMap.initDir(d,:)); hold on;
            tr = trialRanges.initDir{d};
            stX = steadyX.initDir(d,:);
        else % RT
            figure(2);
            subplot(4,8, subplotMap.RT(d,:)); hold on;
            tr = trialRanges.RT{d};
            stX = steadyX.RT(d,:);
        end

        % --- Per-subject traces ---
        for subj = 1:numSubjects
            plot_metric(totSW{subj}.(directions{d}).(metrics{m}), tr, [], [], directions{d}, [], [], colors, col, trialAfterSW, 1, []);
        end

        for subj = 1:numSubjects
             subRTave.AtoB(subj,:) = nanmean(totSW{subj}.AtoB.RT);
             subRTave.BtoA(subj,:) = nanmean(totSW{subj}.BtoA.RT);
        end
                   
        % --- Group average ---
        plot_metric(dCombined.(metrics{m}).(directions{d}), tr, stX, metrics{m}, directions{d}, SS, eMix, colors, col, trialAfterSW, 2, []);

        % Axes and labels
        axis_label(m,metrics,xRanges,yRanges,xtickVals)

    end
end

%% ===== Plot Individuals =====

if any(~cellfun('isempty', {reps.subj}))

    disp('Representative subjects:');
    disp([reps.subj])

    for rsub = 1:numel(reps)

        subj = reps(rsub).subj;

        for m = 1:2
            for d = 1:2 % direction

                if m == 1 % initDir
                    figure(1); subplot(4,8,subplotPos.initDir{rsub}); hold on;
                    tr = trialRanges.initDir{d};
                else % RT
                    figure(2);
                    if d == 1
                        subplot(4,8,subplotPos.RT_BtoA{rsub}); hold on;
                    else
                        subplot(4,8,subplotPos.RT_AtoB{rsub}); hold on;
                    end
                    tr = trialRanges.RT{d};
                end

                plot_metric(totSW{subj}.(directions{d}).(metrics{m}), tr, [], metrics{m}, directions{d}, SS, eMix, colors, col, trialAfterSW, 3, subj);
                axis_label(m,metrics,xRanges,yRanges,xtickVals)

            end
        end
    end
else
    disp('No representative subjects specified. Skipped.');
end

%% ===== Plot policy =====
plot_policy(Result, PR, directions, metrics, reps, trialRanges, subplotMap,subplotPos, xRanges, xtickVals)


%% ===== Functions ===== 

% ---------- Load data ---------- 
% === Load data from Exp1 or Exp2 ===
function [cfg, SW, SS, reps, numSubjects, dataCO, PR] = load_experiment(expID)
% define params(file name, representative subjects, figure number offset, tags for saving) for each experiment
    switch expID
        case 1  % ===== Exp1 =====
            cfg.tag           = 'Exp1';
            cfg.file.switch   = 'Exp1_Switch.mat';
            cfg.file.ss       = 'Exp1_SteadyState.mat';
            cfg.file.polrec   = 'Exp1_PolicyRecovery.mat';
            cfg.figBase       = 100;    % Figure number offset
            reps = struct('subj', {7,11});
        case 2  % ===== Exp2 =====
            cfg.tag           = 'Exp2';
            cfg.file.switch   = 'Exp2_Switch.mat';
            cfg.file.ss       = 'Exp2_SteadyState.mat';
            cfg.file.polrec   = 'Exp2_PolicyRecovery.mat';
            cfg.figBase       = 200;
            reps = struct('subj', {});
        otherwise
            error('Unknown expID: %d', expID);
    end

    % --- Load switching test ---
    S = safe_load(cfg.file.switch);
    SW = pick_field(S, {'_SW','SW','Switch'});   % pick data even variables have slightly different names

    % --- Load steady state ---
    T = safe_load(cfg.file.ss);
    SS = pick_field(T, {'_SS','SS','Steady'});
    dataCO = pick_field(T, {'dco','dco_'});

    % --- Total number of subjects ---
    numSubjects  = size(SW, 2); 

    % --- Policy Recovery --- 
    P = safe_load(cfg.file.polrec);
    PR = pick_field(P, {'Result'});

end

function S = safe_load(fname)
% check if a file exits or not, and then load data to S
    if exist(fname,'file') ~= 2
        error('File not found: %s', fname);
    end
    S = load(fname);
end

function val = pick_field(S, keys)
% get back the first field which has one of keys from S
    fn = fieldnames(S);
    idx = [];
    for i = 1:numel(keys)
        hit = find(contains(fn, keys{i}, 'IgnoreCase',true), 1, 'first');
        if ~isempty(hit), idx = hit; break; end
    end
    if isempty(idx), idx = 1; end           % if it's empty, back the first one
    val = S.(fn{idx});
end
% ----------------------------


% ---------- Policy ---------- 
% === Get distribution from center-out task ===
function [targAng, initDir] = load_steady_maps(sub, dataCO)
    for sw = 1:2
        if sw == 1 % Baseline to De novo
            % Load steady-state De novo
            targAng{sw} = dataCO{sub,2}.targAng; % targetAngle
            initDir{sw,1} = dataCO{sub,2}.initDirL; % hd = 1:2 % 1:LH, 2:RH
            initDir{sw,2} = dataCO{sub,2}.initDirR;
        elseif sw == 2 % De novo to Baseline
            % Load steady-state Baseline
            targAng{sw} = dataCO{sub,1}.targAng; % targetAngle
            initDir{sw,1} = dataCO{sub,1}.initDirL; % hd = 1:2 % 1:LH, 2:RH
            initDir{sw,2} = dataCO{sub,1}.initDirR;
        end
    end
end

% === Edged wraped up ===
% add points to the edges and use different bin size for x and y
% -200 to 200 deg
function [targAng, initDir] = wrap_angles(targAng, initDir)
    for sw = 1:2
        idx_low  = targAng{sw} > -180 & targAng{sw} < -160;
        idx_high = targAng{sw} > 160  & targAng{sw} < 180;
    
        % wrap-around target angles
        extraTarg = [targAng{sw}(idx_low) + 360, targAng{sw}(idx_high) - 360];
    
        % wrap-around initDir for both hands
        for hd = 1:2
            extraInit{hd} = [initDir{sw,hd}(idx_low), initDir{sw,hd}(idx_high)];
            initDir{sw,hd} = [initDir{sw,hd}, extraInit{hd}];
        end
    
        targAng{sw} = [targAng{sw}, extraTarg];
    end
end

% === Build probability map ===
function p_map = build_probability_maps(targAng, initDir, mw, binX, binY)

    % --- Smooting (get moving average) ---
    for sw = 1:2
        for hd = 1:2 % hand, 1:left hand, 2:right hand
            HtargMat{sw,hd} = zeros(100,100);
        end
        for hd = 1:2
            for i = 1:size(targAng{sw},2)
                cntT = 1;
                while -200 + (cntT-1)*360/mw + 360/binX <= 200
                    if  targAng{sw}(i)>-200+(cntT-1)*360/mw && targAng{sw}(i)<=-200+(cntT-1)*360/mw+360/binX && ~isnan(targAng{sw}(i))
                        cntH = 1;
                        while -200 + (cntH-1)*360/mw + 360/binY <= 200
                            if initDir{sw,hd}(i)>-200+(cntH-1)*360/mw && initDir{sw,hd}(i)<=-200+(cntH-1)*360/mw+360/binY && ~isnan(initDir{sw,hd}(i))
                                HtargMat{sw,hd}(cntH,cntT) = HtargMat{sw,hd}(cntH,cntT) + 1;
                            end
                            cntH = cntH + 1;
                        end
                    end
                    cntT = cntT + 1;
                end
            end
        end
    end
    cntT = cntT - 1;
    cntH = cntH - 1;
    
    % --- Normalization ---
    for sw = 1:2
        for hd = 1:2
            colSum = sum(HtargMat{sw,hd}, 1);
            ProHtargMat{sw,hd} = HtargMat{sw,hd} ./ (colSum + (colSum==0));  % avoid div/0
        end
    end
    
    % --- Final policy maps ---
    p_map{1,1} = ProHtargMat{1,1}(1:cntH,1:cntT); % mapA, left
    p_map{1,2} = ProHtargMat{1,2}(1:cntH,1:cntT); % mapA, right
    p_map{2,1} = ProHtargMat{2,1}(1:cntH,1:cntT); % mapB, left
    p_map{2,2} = ProHtargMat{2,2}(1:cntH,1:cntT); % mapB, right

end

% === Estimate maximum liklihood, optimize pA, pB ===
function pAB_opt = estimate_trial_policies(SW, sub, p_map, numtr)

    for sw = 1:2 %1:mapB to mapA, 2:mapA to mapB
    
        % --- theta: all trials post switch ---
        for tr = 1:numtr

            % L:1, R:2, T:3
            % HtargMat{sw,hd}
            if sw == 1
                theta = [];
                theta(:,1) = SW{sub}.BtoA.initDirL(:,tr);
                theta(:,2) = SW{sub}.BtoA.initDirR(:,tr);
                theta(:,3) = SW{sub}.BtoA.targAng(:,tr);
            elseif sw == 2
                theta = [];
                theta(:,1) = SW{sub}.AtoB.initDirL(:,tr);
                theta(:,2) = SW{sub}.AtoB.initDirR(:,tr);
                theta(:,3) = SW{sub}.AtoB.targAng(:,tr);
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
                p_I(n,1) = p_map{1,1}(I(n,1),I(n,3)) * p_map{1,2}(I(n,2),I(n,3)) - (1/CNTT)^2;
                p_I(n,2) = p_map{2,1}(I(n,1),I(n,3)) * p_map{2,2}(I(n,2),I(n,3)) - (1/CNTT)^2;
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

% === Plot policy (group, representative subjects) ===
function plot_policy(Result, PR, directions, metrics, reps, trialRanges, subplotMap, subplotPos, xRanges, xtickVals)
    figure(3); hold on;
    m = 1;

    % Policy recovery
    DE_DE = mean(mean(PR{1,1},2));
    DE_BA = mean(mean(PR{1,2},2));
    DE_N = mean(mean(PR{1,3},2));
    BA_DE = mean(mean(PR{2,1},2));
    BA_BA = mean(mean(PR{2,2},2));
    BA_N = mean(mean(PR{2,3},2));

    % Group
    for d = 1:length(directions) % BaseToDenovo / DenovoToBase
        subplot(4,8, subplotMap.initDir(d,:)); hold on;
        tr = trialRanges.initDir{d};
    
        % Plot probability
        shadedErrorBar(tr,nanmean(Result{d,3}),seNaN(Result{d,3}),{'g-o','markerfacecolor','g','markeredgecolor','g','markersize',6,'LineWidth',2})
        shadedErrorBar(tr,nanmean(Result{d,1}),seNaN(Result{d,1}),{'m-o','markerfacecolor','m','markeredgecolor','m','markersize',6,'LineWidth',2})
        shadedErrorBar(tr,nanmean(Result{d,2}),seNaN(Result{d,2}),{'b-o','markerfacecolor','b','markeredgecolor','b','markersize',6,'LineWidth',2})
    
        % Plot policy recovery
        if d == 1
            plotPolicyRecovery(tr, DE_DE, DE_BA, DE_N)
        else
            plotPolicyRecovery(tr, BA_DE, BA_BA, BA_N);
        end
    end
    policy_axes_labels(xRanges, metrics, m, xtickVals)

    % Representative subjects
    if any(~cellfun('isempty', {reps.subj})) % if there is a representative subject

        disp('Policy representative subjects:');
        disp([reps.subj])

        for rsub = 1:numel(reps)

            subj = reps(rsub).subj;

            for d = 1:2 % direction
                figure(3); subplot(4,8,subplotPos.initDir{rsub}); hold on;
                tr = trialRanges.initDir{d};

                % Plot probability
                plot(tr,Result{d,3}(subj,:),'o-','markerfacecolor','g','markeredgecolor','g','markersize',6,'Color','g','LineWidth',2)
                plot(tr,Result{d,1}(subj,:),'o-','markerfacecolor','m','markeredgecolor','m','markersize',6,'Color','m','LineWidth',2)
                plot(tr,Result{d,2}(subj,:),'o-','markerfacecolor','b','markeredgecolor','b','markersize',6,'Color','b','LineWidth',2)
            end

            policy_axes_labels(xRanges, metrics, m, xtickVals)

        end
    else
        disp('No representative subjects specified. Skipped.');
    end

end

% === Plot policy recovery ===
function plotPolicyRecovery(rngVals, valDE, valBA, valOther)
    plot(rngVals, ones(size(rngVals))*valDE,    'm--','LineWidth',1); % de novo
    plot(rngVals, ones(size(rngVals))*valBA,    'b--','LineWidth',1); % base
    plot(rngVals, ones(size(rngVals))*valOther, 'g--','LineWidth',1); % other
end

% === Plot axes ===
function policy_axes_labels(xRanges, metrics, m, xtickVals)
    axis([xRanges.(metrics{m})(1) xRanges.(metrics{m})(2) 0 1]);
    set(gca,'FontSize',8);
    xlabel('Trial post switch','FontSize',12);
    ylabel('Probability','FontSize',12);
    
    % Switch trial marker
    plot([2 2], [0 1], 'k:');
    
    % X-tick labels
    xticks(0:xRanges.(metrics{m})(2));
    vals = xtickVals.(metrics{m});  % ex: [0 0 0 1 0 0 0 5 ... 10]
    labels = arrayfun(@num2str, vals, 'uni', false);
    labels(vals==0) = {''};  % replace 0 as empty 
    xticklabels(labels);
    xtickangle(0);

end
% ----------------------------


% ---------- Initial direction error and RT ---------- 
% === plot individual data and average ===
function plot_metric(data, tr, stX, metName, dirName, SS, eMix, colors, col, trialAfterSW, Index, subj)

    if Index == 3 % plot representative
        plot(tr, data(:,1:trialAfterSW), 'o', 'MarkerSize', 3, ...
            'MarkerEdgeColor', col, 'MarkerFaceColor', col);
        plot(tr, nanmean(data(:,1:trialAfterSW)), ...
            [colors.(dirName) 'o-'], 'MarkerSize', 3, ...
            'MarkerEdgeColor', colors.(dirName), 'MarkerFaceColor', colors.(dirName), ...
            'LineWidth', 2);
        
        if strcmp(dirName, 'AtoB') % DE to BA
            if strcmp(metName, 'initDir') % initDir
                plot(tr, repmat(SS.(metName)(subj,1),1,length(tr)),'b-', 'LineWidth', 1); % BA
                plot(tr, repmat(SS.initDir_BA(subj,2),1,length(tr)),'m:', 'LineWidth', 1.2); % DE reconstruct under BA
                plot(tr, eMix{2}(subj,:), 'bo:','LineWidth',1.2);
            else % RT
                plot(tr, repmat(SS.(metName)(subj,1),1,length(tr)),'k--'); % BA
            end
        elseif strcmp(dirName, 'BtoA') % BA to DE
            if strcmp(metName, 'initDir') % initDir
                plot(tr, repmat(SS.(metName)(subj,2),1,length(tr)),'m-', 'LineWidth', 1); % DE
                plot(tr, repmat(SS.initDir_DE(subj,1),1,length(tr)),'b:', 'LineWidth', 1.2); % BA reconstruct under DE
                plot(tr, eMix{1}(subj,:), 'mo:','LineWidth',1.2);
            else % RT
                plot(tr, repmat(SS.(metName)(subj,2),1,length(tr)),'k--'); % DE
            end
        end
    
    elseif Index == 1 % Per-subject traces        
        plot(tr,nanmean(data(:,1:trialAfterSW)), '-', 'Color', col)
    
    elseif Index == 2 % Group average
        plot(tr, nanmean(data(:,1:trialAfterSW)), ...
            [colors.(dirName) 'o-'], 'MarkerSize', 3, ...
            'MarkerEdgeColor', colors.(dirName), 'MarkerFaceColor', colors.(dirName),...
            'LineWidth', 2);    

        % Steady state line
        if strcmp(metName, 'initDir') && strcmp(dirName, 'AtoB')
            plot(stX, repmat(nanmean(SS.(metName)(:,1)), 1, 2), 'b-', 'LineWidth', 1); % Steady state BA
            plot(stX, repmat(nanmean(SS.initDir_BA(:,2)), 1, 2), 'm:', 'LineWidth', 1.2); % DE reconstruct under BA
            plot(tr, nanmean(eMix{2}), 'bo:','LineWidth',1.2); 
        elseif strcmp(metName, 'initDir') && strcmp(dirName, 'BtoA') % Base -> Denovo
            plot(stX, repmat(nanmean(SS.(metName)(:,2)), 1, 2), 'm-', 'LineWidth', 1); % Steady state line DE
            plot(stX, repmat(nanmean(SS.initDir_DE(:,1)), 1, 2), 'b:', 'LineWidth', 1.2); % BA reconstruct under DE
            plot(tr, nanmean(eMix{1}), 'mo:','LineWidth',1.2);
        elseif strcmp(metName, 'RT') && strcmp(dirName, 'AtoB')
            plot(tr, repmat(nanmean(SS.(metName)(:,1)), 1, length(tr)),'k--'); % DE
        elseif strcmp(metName, 'RT') && strcmp(dirName, 'BtoA')
            plot(tr, repmat(nanmean(SS.(metName)(:,2)), 1, length(tr)),'k--'); % DE
        end

    end
end


function axis_label(m,metrics,xRanges,yRanges,xtickVals)

    % Axes and labels
    axis([xRanges.(metrics{m})(1) xRanges.(metrics{m})(2) yRanges.(metrics{m})]);
    set(gca,'FontSize',8);
    xlabel('Trial post switch','FontSize',12);
    ylabel(strrep(metrics{m},'_',' '),'FontSize',12);
    
    % Switch trial marker
    plot([2 2], yRanges.(metrics{m}), 'k:');
    
    % X-tick labels
    xticks(0:xRanges.(metrics{m})(2));
    xtickangle(0);

    vals = xtickVals.(metrics{m}); % ex: [0 0 0 1 0 0 0 5 ... 10]
    labels = arrayfun(@num2str, vals, 'uni', false);
    labels(vals==0) = {''}; % if 0, its label is empty
    xticklabels(labels);
 
    if strcmp(metrics{m},'RT')
        yticks(200:200:800);
    end

end
% ----------------------------



%% Figure 3I (Experiment 1)
% Predicted initial direction error vs Absolute initial direction error 

if expID == 1
    nDir = numel(directions);

    % Set a panel
    nCol = ceil(sqrt(trialAfterSW));
    nRow = ceil(trialAfterSW / nCol);

    % Width and height per panle (inch)
    panelWidth  = 3; 
    panelHeight = 3;

    % Get size of whole figure 
    figWidth  = panelWidth * nCol;
    figHeight = panelHeight * nRow;

    % Create Figure（backgroud color: white, size: fixed）
    figure(4); clf; set(gcf,'Color', 'w', 'Units', 'inches', 'Position', [xFigPos yFigPos figWidth figHeight]);

    tiledlayout(nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');

    statsX = []; statsY = [];
    for tr = 1:trialAfterSW
        nexttile; hold on;

        for subj = 1:numSubjects
            for d = 1:numel(directions)
                % x axis: mean initDir per trials post switch
                x = nanmean(totSW{subj}.(directions{d}).initDir(:,tr));
                % y axis: eMix for at that trials
                y = eMix{d}(subj, tr);

                % for stats outside of matlab (get Bayes facotor) 
                if tr == 1
                    statsX(end+1,1) = x;
                    statsY(end+1,1) = y;
                end

                plot(x, y, 'o', ...
                    'MarkerEdgeColor', colors.(directions{d}), ...
                    'MarkerFaceColor', colors.(directions{d}), ...
                    'MarkerSize',4);
            end
        end
        plot([0 80],[0 80], 'k-')

        title(sprintf('Trial %d', tr));
        xlabel('initDir (deg)');
        ylabel('eMix');
        xticks(0:20:80);
        yticks(0:20:80);
        axis equal
        xlim([0 80]);
        ylim([0 80]);
        % grid on;
    end

    % Legend
    lgd = legend(directions, 'Location','bestoutside');

end

%% Figure 2DE, correlation between initial direction error and RT-SS (RT differences over steady state)
% Experiment 1

if expID == 1

    % === Parameters ===
    endtr    = 6;       % post switch tr > 6 -> steady state
    cond     = 1;       % SW{cond,subj}

    % ============= Preprocessing for group average（mRT / minitDir） =============
    % get mRT.(dir)(subj,tr), minitDir.(dir)(subj,tr)
    [mRT, minitDir] = compute_means_per_trial(totSW, cond, numSubjects, trialAfterSW);

    % ============= Figure: 1st trial, 6-10 trials post switch (B->A / A->B) =============
    figure(5); clf; set(gcf,'units','inch','position',[xFigPos,yFigPos,xFigSize,yFigSize/4],'Color','w');

    if ~exist('endtr','var') || isempty(endtr), endtr = 6; end

    % ===== Baseline to De novo (B->A) =====
    subplot(1,2,1); hold on;
    for subj = 1:numSubjects
        % 1st trial post switch
        plot(minitDir.BtoA(subj,1), ...
            mRT.BtoA(subj,1) - SS.RT(subj,2), ...
            'o','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);

        % average of 6–10 trials
        plot(nanmean(minitDir.BtoA(subj,endtr:trialAfterSW)), ...
            nanmean(mRT.BtoA(subj,endtr:trialAfterSW) - SS.RT(subj,2)), ...
            'o','MarkerEdgeColor','m','MarkerSize',5);

        % line: 1st trial ↔ average of 6–10 trials
        plot([minitDir.BtoA(subj,1), nanmean(minitDir.BtoA(subj,endtr:trialAfterSW))], ...
            [mRT.BtoA(subj,1)-SS.RT(subj,2), nanmean(mRT.BtoA(subj,endtr:trialAfterSW)-SS.RT(subj,2))], 'm-');
    end
    % Correlation (1trial only, ΔRT = RT - SS）
    xBA = minitDir.BtoA(:,1);
    yBA = mRT.BtoA(:,1) - SS.RT(:,2);
    % Correlation (6-10 trials, ΔRT = RT - SS）
    % xBA = nanmean(minitDir.BtoA(:,endtr:trialAfterSW),2);
    % yBA = nanmean(mRT.BtoA(:,endtr:trialAfterSW) - SS.RT(:,2),2);
    title(report_corr(xBA, yBA, '1st tr'));
    xlim([0 80]); ylim([-200 300]);
    ylabel('RT - Steady State (ms)','FontSize',12);
    xlabel('Absolute Init Dir Error (deg)','FontSize',12);
    box on;


    % ===== De novo to Baseline (A->B) =====
    subplot(1,2,2); hold on;
    for subj = 1:numSubjects
        % 1st trial post switch
        plot(minitDir.AtoB(subj,1), ...
            mRT.AtoB(subj,1) - SS.RT(subj,1), ...
            'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);

        % average of 6–10 trials
        plot(nanmean(minitDir.AtoB(subj,endtr:trialAfterSW)), ...
            nanmean(mRT.AtoB(subj,endtr:trialAfterSW) - SS.RT(subj,1)), ...
            'o','MarkerEdgeColor','b','MarkerSize',5);

        % line: 1st trial ↔ average of 6–10 trials
        plot([minitDir.AtoB(subj,1), nanmean(minitDir.AtoB(subj,endtr:trialAfterSW))], ...
            [mRT.AtoB(subj,1)-SS.RT(subj,1), nanmean(mRT.AtoB(subj,endtr:trialAfterSW)-SS.RT(subj,1))], ...
            'b-');
    end
    % Correlation (1trial only, ΔRT = RT - SS）
    xAB = minitDir.AtoB(:,1);
    yAB = mRT.AtoB(:,1) - SS.RT(:,1);
    % Correlation (6-10 trials, ΔRT = RT - SS）
    % xAB = nanmean(minitDir.AtoB(:,endtr:trialAfterSW),2);
    % yAB = nanmean(mRT.AtoB(:,endtr:trialAfterSW) - SS.RT(:,1),2);
    title(report_corr(xAB, yAB, '1st tr'));
    xlim([0 80]); ylim([-200 300]);
    ylabel('RT - Steady State (ms)','FontSize',12);
    xlabel('Absolute Init Dir Error (deg)','FontSize',12);

end

% ===== Functions =====
function [mRT, minitDir] = compute_means_per_trial(SW, cond, numSubjects, nTrials)
    directions = {'BtoA','AtoB'};
    for k = 1:numel(directions)
        d = directions{k};
        mRT.(d)      = nan(numSubjects, nTrials);
        minitDir.(d) = nan(numSubjects, nTrials);
        for subj = 1:numSubjects
            D = SW{cond,subj}.(d);
            % average of columns (= average of trials) (ignore NaN in row direction）
            mRT.(d)(subj,:)      = nanmean(D.RT(:,1:nTrials),      1);
            minitDir.(d)(subj,:) = nanmean(D.initDir(:,1:nTrials), 1);
        end
    end
end

function s = report_corr(x, y, label)
    % x,y: vectors with NaNs allowed
    idx = ~(isnan(x) | isnan(y));
    x = x(idx); y = y(idx);
    n = numel(x); df = n - 2;

    [R,P] = corrcoef(x,y); r = R(1,2); p = P(1,2);

    % Fisher z CI
    z  = atanh(r);
    se = 1/sqrt(n-3);
    z95 = 1.96 * se;
    ci_r = tanh([z - z95, z + z95]);

    s = sprintf('%s: r(%d) = %.2f, 95%% CI [%.2f, %.2f], p = %.2f, n = %d', ...
                label, df, r, ci_r(1), ci_r(2), p, n);
end


% return;

%% Exp2, plot participant's choice of mapping (shown by button press)

if expID == 2

    %%%%%%%%%%%%%%%%%%%% draw maps %%%%%%%%%%%%%%%%%%%%
    % Define bins
    edgesT = -185:10:195;   % → bins: [-185,-175), ..., [185,195)
    edgesH = -170:10:200;   % → bins: [-170,-160), ..., [190,200)
    bintarg = 39;

    % load p_map_all.mat
    load Exp2_Switch.mat

    % subjects: who will be analyzed? (can edit depending on purporses)
    if exist('DExp2','var') && iscell(DExp2)
        subList = 1:numel(DExp2);   % All subjects
    else
        error('Cannot find DExp2.');
    end

    nSub = numel(subList);

    res.no_ino_bno = zeros(nSub); res.no_ino_bsw = zeros(nSub); res.no_isw_bno = zeros(nSub); res.no_isw_bsw = zeros(nSub);
    res.sw_ino_bno = zeros(nSub); res.sw_ino_bsw = zeros(nSub); res.sw_isw_bno = zeros(nSub); res.sw_isw_bsw = zeros(nSub);

    % focus on switch/switch trials
    res.BtoD_base = zeros(nSub);  res.BtoD_denovo = zeros(nSub);  res.BtoD_other = zeros(nSub);
    res.DtoB_base = zeros(nSub);  res.DtoB_denovo = zeros(nSub);  res.DtoB_other = zeros(nSub);

    for si = 1:nSub

        sub = subList(si);
        Dsub = DExp2{sub};

        % ---- Angle -> bin (vectorization) ----
        thetaL = Dsub.initDirL(:);
        thetaR = Dsub.initDirR(:);
        thetaT = Dsub.targAng(:);

        I3 = discretize(thetaT, edgesT);   % Target
        I1 = discretize(thetaL, edgesH);   % Left hand
        I2 = discretize(thetaR, edgesH);   % Right hand

        % Smaller than the minimum -> bin 1
        I3(isnan(I3) & thetaT < edgesT(1)) = 1;
        I1(isnan(I1) & thetaL < edgesH(1)) = 1;
        I2(isnan(I2) & thetaR < edgesH(1)) = 1;

        % Bigger than the maximum -> NaN
        I3(thetaT >= edgesT(end)) = NaN;   % Exclude over 195
        I1(thetaL >= edgesH(end)) = NaN;   % Exclude over 200
        I2(thetaR >= edgesH(end)) = NaN;

        for tr = 1:length(I1)
            if ~isnan(I1(tr)) && ~isnan(I2(tr)) && ~isnan(I3(tr))
                p_I(tr,1) = p_map_all{sub}{1,1}(I1(tr),I3(tr))*p_map_all{sub}{1,2}(I2(tr),I3(tr)); % de novo
                p_I(tr,2) = p_map_all{sub}{2,1}(I1(tr),I3(tr))*p_map_all{sub}{2,2}(I2(tr),I3(tr)); % base
                p_I(tr,3) = (1/bintarg) * (1/bintarg);
            else
                p_I(tr,:) = NaN;
            end

            if p_I(tr,1) > p_I(tr,2) && p_I(tr,1) > p_I(tr,3)
                p_I(tr,4) = 1; % de novo
            elseif p_I(tr,2) > p_I(tr,1) && p_I(tr,2) > p_I(tr,3)
                p_I(tr,4) = 2; % base
            else
                p_I(tr,4) = 100;
            end

        end

        %%%%%% plot true event/participants' choice %%%%%%
        % divide data into block
        block = find(DExp2{sub}.tFile(:,1)==1);
        for h = 1:length(block)

            if h < length(block)
                strp = block(h);
                endp = block(h+1)-1;
            elseif h == length(block)
                strp = block(h);
                endp = length(DExp2{sub}.tFile(:,1));
            end

            D{h}.tFile = DExp2{sub}.tFile(strp:endp,:);
            D{h}.p_I = p_I(strp:endp,:);

        end

        swtr = [];
        for blk = 1:size(block,1)
            temp = []; temp2 = [];
            temp = [0; diff(D{blk}.tFile(:,6))]; % check true event changes, add 0 for the 1st trial

            % true event: diff(D{blk}.tFile(:,6))=-1/1, 0(switch, no switch)
            % bottun press (subject's intention): D{blk}.tFile(:,7)=1 (1:press, 0:no)
            % chosen policy (subject's behaviour): D{blk}.p_I(:,4)=1,2,100(de novo, base, other)
            % true map: D{blk}.tFile(:,6)=0,1(base, de novo)
            temp2 = [temp, D{blk}.tFile(:,7), D{blk}.p_I(:,4), D{blk}.tFile(:,6)];

            % add all blocks
            swtr = [swtr; temp2];
        end

        % true event / intention / behavior
        for i = 2:length(swtr(:,1))
            if swtr(i,1) == 0 && swtr(i,2) == 0 && swtr(i,3) == 0     % no switch / no switch / no switch
                res.no_ino_bno(sub) = res.no_ino_bno(sub) + 1;
            elseif swtr(i,1) == 0 && swtr(i,2) == 0 && swtr(i,3) ~= 0 % no switch / no switch / switch
                res.no_ino_bsw(sub) = res.no_ino_bsw(sub) + 1;
            elseif swtr(i,1) == 0 && swtr(i,2) == 1 && swtr(i,3) == 0 % no switch / switch / no switch
                res.no_isw_bno(sub) = res.no_isw_bno(sub) + 1;
            elseif swtr(i,1) == 0 && swtr(i,2) == 1 && swtr(i,3) ~= 0 % no switch / switch / switch
                res.no_isw_bsw(sub) = res.no_isw_bsw(sub) + 1;
            elseif swtr(i,1) ~= 0 && swtr(i,2) == 0 && swtr(i,3) == 0 % switch / no switch / no switch
                res.sw_ino_bno(sub) = res.sw_ino_bno(sub) + 1;
            elseif swtr(i,1) ~= 0 && swtr(i,2) == 0 && swtr(i,3) ~= 0 % switch / no switch / switch
                res.sw_ino_bsw(sub) = res.sw_ino_bsw(sub) + 1;
            elseif swtr(i,1) ~= 0 && swtr(i,2) == 1 && swtr(i,3) == 0 % switch / switch / no switch
                res.sw_isw_bno(sub) = res.sw_isw_bno(sub) + 1;
            elseif swtr(i,1) ~= 0 && swtr(i,2) == 1 && swtr(i,3) ~= 0 % switch / switch / switch
                res.sw_isw_bsw(sub) = res.sw_isw_bsw(sub) + 1;
            end
        end
        res.SW(sub) = res.sw_ino_bno(sub) + res.sw_ino_bsw(sub) + res.sw_isw_bno(sub) + res.sw_isw_bsw(sub);
        res.NO(sub) = res.no_ino_bno(sub) + res.no_ino_bsw(sub) + res.no_isw_bno(sub) + res.no_isw_bsw(sub);

        % Focus on bahaviour at Switch(true event)/Switch(intention) trials
        for i = 2:length(swtr(:,1))
            if swtr(i,1) == 1 && swtr(i,2) == 1 % true: Base -> De novo / intention: switch
                if swtr(i,3) == 1 % behavior: de novo
                    res.BtoD_denovo(sub) = res.BtoD_denovo(sub) + 1;
                elseif swtr(i,3) == 2 % behavior: base
                    res.BtoD_base(sub) = res.BtoD_base(sub) + 1;
                elseif swtr(i,3) == 100 % behavior: something else
                    res.BtoD_other(sub) = res.BtoD_other(sub) + 1;
                end
            elseif swtr(i,1) == -1 && swtr(i,2) == 1 % true: De novo -> Base / intention: switch
                if swtr(i,3) == 1 % behavior: de novo
                    res.DtoB_denovo(sub) = res.DtoB_denovo(sub) + 1;
                elseif swtr(i,3) == 2 % behavior: de novo
                    res.DtoB_base(sub) = res.DtoB_base(sub) + 1;
                elseif swtr(i,3) == 100 % behavior: something else
                    res.DtoB_other(sub) = res.DtoB_other(sub) + 1;
                end
            end
        end
        res.BtoD(sub) = res.BtoD_base(sub) + res.BtoD_denovo(sub) + res.BtoD_other(sub);
        res.DtoB(sub) = res.DtoB_base(sub) + res.DtoB_denovo(sub) + res.DtoB_other(sub);

    end

    xFigPos = 1200; yFigPos = 200; figWidth = 800; figHeight = 300;
    figure(100); clf; set(gcf,'Color', 'w', 'Units', 'points', 'Position', [xFigPos yFigPos figWidth figHeight]);
    hold on;

    % group average: true event/participants' choice
    allSW = sum(res.SW(:));
    allNO = sum(res.NO(:));
    subplot(1,2,1); hold on;
    x = ["No/No" "No/SW" "Sw/No" "Sw/Sw"];
    y = [(sum(res.no_ino_bno(:))+sum(res.no_ino_bsw(:)))/allNO, (sum(res.no_isw_bno(:))+sum(res.no_isw_bsw(:)))/allNO, ...
        (sum(res.sw_ino_bno(:))+sum(res.sw_ino_bsw(:)))/allSW, (sum(res.sw_isw_bno(:))+sum(res.sw_isw_bsw(:)))/allSW];
    bar(x,y,0.4)
    set(gca,'FontSize',16);
    ylabel('%','FontSize',16)
    ylim([0 1])
    title('True Event / Intention','Fontsize',16)

    % group average: mapping chosen at switch trials (switch/switch)
    allDtoB = sum(res.DtoB(:));
    allBtoD = sum(res.BtoD(:));
    subplot(1,2,2); hold on;
    x = ["Base to De Novo" "De Novo to Base"];
    y = [sum(res.BtoD_base(:))/allBtoD sum(res.BtoD_denovo(:))/allBtoD sum(res.BtoD_other(:))/allBtoD;
        sum(res.DtoB_base(:))/allDtoB sum(res.DtoB_denovo(:))/allDtoB sum(res.DtoB_other(:))/allDtoB];
    b = bar(x,y);
    b(1).FaceColor = 'b'; % base
    b(2).FaceColor = 'm'; % de novo
    b(3).FaceColor = 'g'; % other

    set(gca,'FontSize',16);
    ylabel('%','FontSize',16)
    ylim([0 1])
    title('Mapping chosen at SW/SW trials','Fontsize',16)

end





