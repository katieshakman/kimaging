% Plot Multiple AvgZscores

clear; 
close all; 
cd('/Volumes/KATIELAB1/Data_Analysis/B-1_Stim/RealMV1-296B-L238_B-1_Pairing')
%% Set labels and user parameters
xLab = 'Seconds'; 
yLab = 'Z-Score'; 
stimOn = 6; 
stimOff = 8; 
totalTime = 15; 
auc_type = 'stimToEnd';

dropFrames =0; % Set to 1 in order to have script look for frames when stim was on and drop those, set to 2 to also drop each frame after each of those in case the light was not fully gone when shutter opened. 
% framesToDrop = [51:2:63]; % 
% framesToDrop = [32:1:40]; % 901 Live 18Bstim e.g. 12/1/2017
% framesToDrop = 32; 
% framesToDrop = 32:2:40; 
% framesToDrop = [52:2:64]; % 901 explants, 7 pulses
% framesToDrop = 52; % 901 explants, 1 pulse
% framesToDrop = [25:2:31]; % 901 Live MB050B and MB543B stim
%% Choose data to load
%files = V1cV1g_loadDataFor_PlotMultipleAvgZscores_ListInput('expt'); textname = 'V1cV1g BitterLight';

group='expt_post_CS-'; 
files = L238MV1cV1gLive_loadDataFor_PLotMultipleAvgZscores_ListInput(group);
textname = 'L238MV1cV1gLive';

% group='ctrl2u_post_CS+'; 
% files = MV1cV1gLive_loadDataFor_PLotMultipleAvgZscores_ListInput(group);
% textname = 'MV1cV1gLive';

% group='ctrl1_pre_CS-'; 
% files = MV1cV1gLive_loadDataFor_PLotMultipleAvgZscores_ListInput(group);
% textname = 'MV1cV1gLive';

% group = 'ctrl1'; 
% files = SS01194cV1g7pulse_loadDataFor_PlotMultipleAvgZscores_ListInput(group); 
% textname = 'a2scSS01194cV1gExplants7pulses'; 

% files = MB018BcV1g_loadDataFor_PlotMultipleAvgZscores_ListInput('expt'); textname = 'MB018BcV1g'; %input to function is 'ctrl1' or 'expt'
% files = MB018BcMV1g_loadDataFor_PlotMultipleAvgZscores_ListInput('expt'); textname = 'MB018BcMV1g'; %input to function is 'ctrl1' or 'expt'
% files = MB018BcV1g_loadDataFor_PlotMultipleAvgZscores_ListInput('ctrl1'); textname = 'MB018BcV1g';

% files = MB018BcV1gLiveRm901_loadFor_PlotMultAvgZscores_ListInput('expt'); textname = 'MB018BcV1gLiveRm901';

% group='ctrl2u'; files =
% MB018BcV1gExplant_loadDataFor_PlotMultipleAvgZscores_ListInput(group);
% textname = 'MB018cV1gExplant'; % group is 'ctrl1' or 'expt' or 'ctrl2u'

% group='expt_pre_CS+'; 
% files = THcV1gLive_loadDataFor_PlotMultipleAvgZscores_ListInput(group);
% textname = 'THcV1gLive'; 

% group='expt_post_CS+'; 
% files = THcV1gLive_loadDataFor_PlotMultipleAvgZscores_ListInput(group);
% textname = 'THcV1gLive';

% group='ctrl1_post_CS+'; 
% files = PAMcV1gLive_loadDataFor_PLotMultipleAvgZscores_ListInput(group);
% textname = 'PAMcV1gLive';

% group='ctrl2u_post_CS-'; 
% files = PAMcV1gLive_loadDataFor_PLotMultipleAvgZscores_ListInput(group);
% textname = 'PAMcV1gLive';

% group='expt_post_CS+'; 
% files = PAMcV1gLive_loadDataFor_PLotMultipleAvgZscores_ListInput(group);
% textname = 'PAMcV1gLive';
 
% group='expt_post_CS-'; 
% files = PAMcV1gLive_loadDataFor_PLotMultipleAvgZscores_ListInput(group);
% textname = 'PAMcV1gLive'; 

% group='ctrl1'; 
% files = MB077CcMV1g7pulse_loadDataFor_PlotMultipleAvgZscores_ListInput(group);
% textname = 'MB077CcMV1g7pulses'; % group is 'ctrl1' or 'expt' or 'ctrl2u'

% group='ctrl1'; 
% files = MB112cV1g_loadDataFor_PlotMultipleAvgZscores_ListInput(group);
% textname = 'MB112CcV1g3pulses'; % group is 'ctrl1' or 'expt' or 'ctrl2u'

% group = 'expt_light_Mch'; 
% files = MB050BcV1gLive_loadDataFor_PLotMultipleAvgZscores_ListInput(group); 
% textname = 'MB050BcV1gLiveDropFrames'; 

% group = 'expt_mch_pre'; 
% files = MB077CcV1gLive7_loadDataFor_PlotMultipleAvgZscores_ListInput(group); 
% textname = 'MB077CcV1gLive'; 

% group = 'ctrl1_post_Oct'; 
% files = MB0543BcV1gLive_loadDataFor_PLotMultipleAvgZscores_ListInput(group);
% textname = 'MB0543BcV1gLive'; 

% group = 'ctrl1'; 
% files = V1cV1g_loadDataFor_PlotMultipleAvgZscores_ListInput(group); 
% textname = 'MB058BcV1gLiveBitter'; 

%% Load data and process data and make first plot
figure; whitebg('w');
timeLen = 0; % initialize
longestTimes = cell(1,1); 
traces = cell(length(files),1); cd
times = cell(length(files),1); 
hold on; 
for filei = 1:length(files)
    filename = files{filei,1};
    load(filename); % loads and timePoints and avgZscores
    disp(strcat('Loaded:',filename)); 
    if length(timePoints)>timeLen
        timeLen = length(timePoints); 
        longestTimes = timePoints;
    end
    traces{filei,1} = avgZscores; 
    times{filei,1} = timePoints; 
    plot(timePoints,avgZscores,'LineWidth',5); 
end

% Make a matrix of the traces, for inspection purposes:
tracesMat = cell2mat(traces); 

% Drop frames if specified by user parameters: 
tracesOrig = tracesMat; % Don't drop any frames from orig reference copy.
if dropFrames == 1 % drop only the specified frames from dFoF
    tracesMat(:,framesToDrop) = nan;
elseif dropFrames == 2
    tracesMat(:,framesToDrop) = nan;
    tracesMat(:,framesToDrop+1) = nan;
elseif dropFrames == 0
    assert(isequaln(tracesMat,tracesOrig)); % do nothing; assert tracesMat still equal to tracesOrig
end

%% Plot the traces for each fly and the average
% Fill in nan values via interpolation: 
% nanx = isnan(x);
% t    = 1:numel(x);
% x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
tracesNew = zeros(size(tracesMat)); 
traceLabels = cell(1,length(files)+1); 
for filei = 1:length(files)
    x = tracesMat(filei,:); 
    nanx = isnan(x); 
    t = 1:numel(x); 
    x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
    tracesNew(filei,:) = x; 
    traceLabels{filei} = ['Fly', num2str(filei)]; 
end
traceLabels{end} = 'Average'; 
% tracesNew = nan(length(traces),length(longestTimes)); 
% for filei = 1:length(files)
%     tracesNew(filei,:) = interp1(times{filei,1},traces{filei,1},longestTimes); 
% end

% Get and plot averages: 
avgZscoreManyFlies = mean(tracesNew,1); 
try
    plot(longestTimes,avgZscoreManyFlies,'LineWidth',5,'Color','k'); 
    xlabel(xLab);  ylabel(yLab);
    title('Trace For Each Fly and Average Across Flies');
    
% lgd = legend('1','2','3','4','5','6','7','Average');
% lgd.Title.String = 'Flies';
catch
    disp('Something went wonky!'); 
end
lgdMany = legend(traceLabels);

%% Make a figure with the mean and SEM error bars for all the traces in this group/folder: 
ntrials = size(tracesNew,1); 
SDdfEachTimepoint = std(tracesNew,1);
SEMdfEachTimepoint = SDdfEachTimepoint./sqrt(ntrials); 
figure; 
try
    plot(longestTimes,avgZscoreManyFlies,'Color','w','LineWidth',3);
    hold on
    errorbar(longestTimes,avgZscoreManyFlies,SEMdfEachTimepoint,'Color','w','LineWidth',3);
    hold off
catch
    disp('Hmm...too few traces to plot with errorbar?'); 
    try
        plot(longestTimes,avgZscoreManyFlies);
%         plot(longestTimes,avgZscoreManyFlies,'Color','w','LineWidth',3);
    catch
        disp('Something went hella wonky!');
    end
end
xlabel(xLab); ylabel(yLab);
title('Average and SEM Across Flies');
legName = [textname, ' n=', num2str(filei)];
lgd = legend(legName);
    
whitebg('k');
ylim([-2 16]); 
xlim([0 ceil(max(longestTimes))]); 


%% Save avg and SEM and max for each fly (after shuttered values removed). 
%  Also save min in the stim interval.
stimOnFr = min(find(longestTimes>stimOn)); 
stimOffFr = min(find(longestTimes>stimOff)); 
minInStimVals = min(tracesMat(:,stimOnFr:stimOffFr),[],2); 

% Get aucs_from_zscore 
aucs_from_zscore = nan(size(tracesMat,1), 1); 
for idx = 1:length(files)
    if strcmp(auc_type,'stimToEnd')
        aucs_from_zscore(idx,1) = trapz(tracesNew(idx,stimOnFr:end)); 
    else
        auc_type = 'stimOnly';
        aucs_from_zscore(idx,1) = trapz(tracesNew(idx,stimOnFr:stimOffFr)); 
    end
end

maxVals = max(tracesMat,[],2);
saveName = [textname,'_',group]; 
save(saveName, 'avgZscoreManyFlies', 'SEMdfEachTimepoint','maxVals','minInStimVals','aucs_from_zscore'); 
display(pwd)


