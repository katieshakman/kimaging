% Plot Multiple AvgZscores

clear; 
close all; 
%% Set labels and user parameters
xLab = 'Seconds'; 
yLab = 'Z-Score'; 

dropFrames = 0; % Set to 1 in order to have script look for frames when stim was on and drop those, set to 2 to also drop each frame after each of those in case the light was not fully gone when shutter opened. 
% framesToDrop = [51:2:63]; % 
% framesToDrop = [32:1:40]; % 901 Live 18Bstim e.g. 12/1/2017
% framesToDrop = 32; 
% framesToDrop = 32:2:40; 
% framesToDrop = [52:2:64]; % 901 explants, 7 pulses
% framesToDrop = 52; % 901 explants, 1 pulse


%% Choose data to load
%files = V1cV1g_loadDataFor_PlotMultipleAvgZscores_ListInput('expt'); textname = 'V1cV1g BitterLight';

% files = MB018BcV1g_loadDataFor_PlotMultipleAvgZscores_ListInput('expt'); textname = 'MB018BcV1g'; %input to function is 'ctrl1' or 'expt'
% files = MB018BcMV1g_loadDataFor_PlotMultipleAvgZscores_ListInput('expt'); textname = 'MB018BcMV1g'; %input to function is 'ctrl1' or 'expt'
% files = MB018BcV1g_loadDataFor_PlotMultipleAvgZscores_ListInput('ctrl1'); textname = 'MB018BcV1g';

% files = MB018BcV1gLiveRm901_loadFor_PlotMultAvgZscores_ListInput('expt'); textname = 'MB018BcV1gLiveRm901';

% group='ctrl2u'; files =
% MB018BcV1gExplant_loadDataFor_PlotMultipleAvgZscores_ListInput(group);
% textname = 'MB018cV1gExplant'; % group is 'ctrl1' or 'expt' or 'ctrl2u'

% group='expt_pre_CS+'; 
% files = THcV1gLive_loadDataFor_PlotMultipleAvgZscores_ListInput(group)
% textname = 'THcV1gLive'; 

group='expt_post_CS+'; 
files = THcV1gLive_loadDataFor_PlotMultipleAvgZscores_ListInput(group);
textname = 'THcV1gLive'; 

% group='ctrl1'; 
% files = MB077CcMV1g7pulse_loadDataFor_PlotMultipleAvgZscores_ListInput(group);
% textname = 'MB077CcMV1g7pulses'; % group is 'ctrl1' or 'expt' or 'ctrl2u'

% group='ctrl1'; 
% files = MB112cV1g_loadDataFor_PlotMultipleAvgZscores_ListInput(group);
% textname = 'MB112CcV1g3pulses'; % group is 'ctrl1' or 'expt' or 'ctrl2u'

%% Load data and process data and make first plot
figure; 
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

%% Fill in nan values via interpolation: 
% nanx = isnan(x);
% t    = 1:numel(x);
% x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
tracesNew = zeros(size(tracesMat)); 
for filei = 1:length(files)
    x = tracesMat(filei,:); 
    nanx = isnan(x); 
    t = 1:numel(x); 
    x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
    tracesNew(filei,:) = x; 
end
% tracesNew = nan(length(traces),length(longestTimes)); 
% for filei = 1:length(files)
%     tracesNew(filei,:) = interp1(times{filei,1},traces{filei,1},longestTimes); 
% end

%% Get and plot averages: 
avgZscoreManyFlies = mean(tracesNew,1); 
try
    plot(longestTimes,avgZscoreManyFlies,'LineWidth',5,'Color','k'); 
xlabel(xLab);  ylabel(yLab); 
title('Trace For Each Fly and Average Across Flies'); 
lgd = legend(textname);
% lgd = legend('1','2','3','4','5','6','7','Average');
% lgd.Title.String = 'Flies';
end
%% Make a figure with the mean and SEM error bars for all the traces in this group/folder: 
ntrials = size(tracesNew,1); 
SDdfEachTimepoint = std(tracesNew,1);
SEMdfEachTimepoint = SDdfEachTimepoint./sqrt(ntrials); 
figure; 
try
    plot(longestTimes,avgZscoreManyFlies,'Color','k','LineWidth',3);
    hold on
    errorbar(longestTimes,avgZscoreManyFlies,SEMdfEachTimepoint,'Color','k','LineWidth',3);
    hold off
    xlabel(xLab); ylabel(yLab);
    title('Average and SEM Across Flies');
    legName = [textname, ' n=', num2str(filei)];
    lgd = legend(legName);
end

whitebg('w');
ylim([-2 5]); 
xlim([2 10]); 

%% Save avg and SEM and max for each fly (after shuttered values removed)
maxVals = max(tracesMat,[],2);
saveName = [textname,'_',group]; 
save(saveName, 'avgZscoreManyFlies', 'SEMdfEachTimepoint','maxVals'); 
display(pwd)

