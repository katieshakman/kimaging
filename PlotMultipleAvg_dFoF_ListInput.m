% Plot Multiple AvgZscores

clear; 
close all; 
% cd('/Volumes/KATIELAB1/Data_Analysis/B-1_Stim/RealMV1-296B-L238_B-1_Pairing')
%% Set labels and user parameters
xLab = 'Seconds'; 
yLab = 'dFoF'; 
stimOn = 6; 
stimOff = 8; 
totalTime = 10; % 15
auc_type = 'stimToEnd';

dropFrames = 1; % Set to 1 in order to have script look for frames when stim was on and drop those, set to 2 to also drop each frame after each of those in case the light was not fully gone when shutter opened. 
% framesToDrop = [51:2:63]; % 
% framesToDrop = [32:1:40]; % 901 Live 18Bstim e.g. 12/1/2017
% framesToDrop = 32; 
% framesToDrop = 32:2:40; 
% framesToDrop = [52:2:64]; % 901 explants, 7 pulses
% framesToDrop = 52; % 901 explants, 1 pulse
framesToDrop = [25:2:31]; % 901 Live MB050B and MB543B stim
%% Choose data to load

% cd('/Users/katherineshakman/CloudStation/Data_Analysis/B-1_Stim/MB027Bstim_B-1_explants_9thFl'); 
cd('/Users/katherineshakman/CloudStation/Data_Analysis/B-1_Stim/MB027Bstim_B-1_explants_9thFl/Data_MB027BcV1gExplant_120pockels_4pulse3vDrop'); 
group='ctrl1';
files = MB027BcV1gExplant_loadDataFor_PlotMultipleAvg_dFoF_ListInput(group); 
% files = MB027BcV1gExpl1pulse1pt3vDr_load_ListInput(group); 
textname = 'MB027BcV1gExplant_dFoF';

%% Load data and process data and make first plot
figure; whitebg('w');
timeLen = 0; % initialize
longestTimes = cell(1,1); 
traces = cell(length(files),1); cd
times = cell(length(files),1); 
hold on; 
for filei = 1:length(files)
    filename = files{filei,1};
    load(filename); % loads avg_dFoF
    traces{filei,1} = avg_dFoF; 
    disp(strcat('Loaded:',filename)); 
    timefilename = strrep(filename,'avg_dFoF','tPts'); 
    try
        load(timefilename);
    catch
        timefilename = strcat('fixed',timefilename); 
    end
    timePoints = tPts(1,:); 
    if length(timePoints)>timeLen
        timeLen = length(timePoints); 
        longestTimes = timePoints;
    end
    times{filei,1} = timePoints; 
    disp(strcat('Loaded times from:',timefilename));
end

% Make a matrix of the traces, for inspection purposes:
smallestLength = Inf;
try 
    tracesMat = cell2mat(traces); 
catch
    % get all the lengths, and take the smallest one
    for idx = 1:length(traces)
        if length(traces{idx}) < smallestLength
            smallestLength = length(traces{idx}); 
        end
    end
    % truncate each trace
    for idx = 1:length(traces)
        tracesShort{idx} = traces{idx}(1:smallestLength); 
    end
    tracesShort = tracesShort'; 
    tracesMat = cell2mat(tracesShort); 
end

% Make initial plot with each trace
for filei = 1:length(files)
    if isinf(smallestLength)
        plot(times{filei,1},traces{filei,1},'LineWidth',5);
    else
        plot(times{filei,1}(1:smallestLength),traces{filei,1}(1:smallestLength),'LineWidth',5);
    end
end
ylabel('dF/F_0'); xlabel('Time (Seconds)'); 

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


%% Add the average trace to the initial plot: 
whitebg('k'); 
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

% Get averages: 
avg_dFoF_ManyFlies = mean(tracesNew,1); 

% Get sem for each fly: 
ntrials = size(tracesNew,1); 
SDdfEachTimepoint = std(tracesNew,1);
SEMdfEachTimepoint_dFoF = SDdfEachTimepoint./sqrt(ntrials); 


% Make sure we are only plotting up to the length of time that the shortest
% trace contains. 
if ~isinf(smallestLength)
   longestTimes = longestTimes(1:smallestLength); 
   avg_dFoF_ManyFlies = avg_dFoF_ManyFlies(1:smallestLength); 
   SEMdfEachTimepoint_dFoF = SEMdfEachTimepoint_dFoF(1:smallestLength); 
end

% Add total average trace to plot: 
try
    if isinf(smallestLength)
        plot(longestTimes,avg_dFoF_ManyFlies,'LineWidth',5,'Color','w');
    else
        plot(longestTimes(1:smallestLength),avg_dFoF_ManyFlies(1:smallestLength),'LineWidth',5,'Color','w');
    end
    xlabel(xLab);  ylabel(yLab);
    title('Trace For Each Fly and Average Across Flies');
    
% lgd = legend('1','2','3','4','5','6','7','Average');
% lgd.Title.String = 'Flies';
catch
    disp('Something went wonky!'); 
end
lgdMany = legend(traceLabels);

%% Make a figure with the mean and SEM error bars for all the traces in this group/folder: 

figure; 
whitebg('w');
try
    plot(longestTimes,avg_dFoF_ManyFlies,'Color','k','LineWidth',3);
    hold on
    errorbar(longestTimes,avg_dFoF_ManyFlies,SEMdfEachTimepoint_dFoF,'Color','k','LineWidth',3);
    hold off
catch
    disp('Hmm...too few traces to plot with errorbar?'); 
    try
        plot(longestTimes,avg_dFoF_ManyFlies);
%         plot(longestTimes,avgZscoreManyFlies,'Color','w','LineWidth',3);
    catch
        disp('Something went hella wonky!');
    end
end
xlabel(xLab); ylabel(yLab);
title('Average and SEM Across Flies');
legName = [textname, ' n=', num2str(filei)];
lgd = legend(legName);

% ylim([-2 16]); 
% xlim([0 ceil(max(longestTimes))]); 

% Add dots to indicate the stimuli
if exist('framesToDrop', 'var')   
    ylimits = ylim; 
    hold on
    xStimMarkers = mean([longestTimes(framesToDrop);longestTimes(framesToDrop+1)],1);
    stimHeights = (-.01)*ones(size(framesToDrop)); 
    scatter(xStimMarkers, stimHeights,40,'filled','k');
end

% % Add line to indicate the stimulus   
% ylimits = ylim;
% hold on
% xStimMarkers = [8,10];
% stimHeights = (-1)*ones(size(xStimMarkers));
% % plot(xStimMarkers, stimHeights,'k','LineWidth',4);
% rectangle('Position', [8,-1.5,2,1],...
%     'FaceColor', .8*[1 1 1], 'EdgeColor','k'); 
% text(8.25,-1,'Bitter On'); 

%% Save avg and SEM and max for each fly (after shuttered values removed). 
%  Also save min in the stim interval.
stimOnFr = min(find(longestTimes>stimOn)); 
stimOffFr = min(find(longestTimes>stimOff)); 
minInStimVals_dFoF = min(tracesMat(:,stimOnFr:stimOffFr),[],2); 

% Get aucs_from_zscore 
aucs_from_dFoF = nan(size(tracesMat,1), 1); 
for idx = 1:length(files)
    if strcmp(auc_type,'stimToEnd')
        aucs_from_dFoF(idx,1) = trapz(tracesNew(idx,stimOnFr:end)); 
    else
        auc_type = 'stimOnly';
        aucs_from_zscore_dFoF(idx,1) = trapz(tracesNew(idx,stimOnFr:stimOffFr)); 
    end
end

maxVals_dFoF = max(tracesMat,[],2);
saveName = [textname,'_',group]; 
save(saveName, 'avg_dFoF_ManyFlies', 'SEMdfEachTimepoint_dFoF','maxVals_dFoF','minInStimVals_dFoF','aucs_from_dFoF'); 
display(pwd)


