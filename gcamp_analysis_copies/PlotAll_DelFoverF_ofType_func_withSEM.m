%% Plot Only Fluor Analysis
% Plot from previously stored deltaF/F data and timepoints. 
% clear all; close all;

%function [] = PlotAll_DelFoverF_ofType_func(imDir)

StartDir = pwd; 

% Variables user can change: 
x_limits = [0,7];
onTime = 3; % time of odor onset
offTime = 5; % time of odor offset

addpath(genpath('/Users/katherineshakman/Documents/MATLAB/raacampbell-shadedErrorBar-1ea286c')); 
addpath(genpath('/Users/katherineshakman/Dropbox/Axel Lab Dropbox/Live MB Imaging Project/B-1_Stim_Backup20180417')); 
addpath(genpath('/Users/katieshak/Documents/MATLAB/raacampbell-shadedErrorBar-1ea286c')); 
addpath(genpath('/Users/katieshak/Dropbox/Axel Lab Dropbox/Live MB Imaging Project/B-1_Stim_Backup20180417')); 
%% Load data files from several user-selected directories. 
% imDir = uigetdir(); % Allow user to select directory. 
% pwdir = pwd; 
% pwdir = pwdir(end-length(imDir)+1:end);
% if ~strcmp(pwdir,imDir)
%     cd(imDir);
% end
imDir = pwd; 
currentDirectory = pwd; % For use in save paths. 

all_delFoverF_files = dir('*delFoverF.mat');
all_delFoverF_files = flipud(all_delFoverF_files);
all_timePoints_files = dir('*timePoints.mat');
all_timePoints_files = flipud(all_timePoints_files);
all_baseline_files = dir('*baselineF0.mat'); 
all_baseline_files = flipud(all_baseline_files);

dFoF = [];
tPts = [];
F0s = []; 
stringsList = []; % initialize
% dFoF = zeros(1, max(size(all_delFoverF_files))); % Initialize
% tPts = zeros(1, max(size(all_timePoints_files))); % Initialize
pattern = '._' ; 
for idx = 1:max(size(all_delFoverF_files))
%     new_delFoverF = load(all_delFoverF_files(idx,1).name);
    try
        notTrueMAT = strfind(all_delFoverF_files(idx,1).name, pattern);
    catch
        notTrueMAT = NaN;
        notTrueMAT = false;
    end
    if notTrueMAT
        ; % do nothing for now.  maybe delete these files in future. 
    else
        load(all_delFoverF_files(idx,1).name);
        dFoF = [delFoverF; dFoF];
        stringsList = [all_delFoverF_files(idx,1).name(9:25); stringsList];
        % used only part of the name that contains date and tseries number
    end
end
stringsList = cellstr(stringsList); % make cell array list of strings
clear idx;
for idx = 1:max(size(all_timePoints_files))
    notTrueMAT = strfind(all_timePoints_files(idx,1).name, pattern);
    if notTrueMAT
        ; % do nothing for now.  maybe delete these files in future. 
    else
        load(all_timePoints_files(idx,1).name);
        tPts = [timePoints; tPts]; 
    end
end
clear idx;
for idx = 1:max(size(all_baseline_files)) % also get baselines
    notTrueMAT = strfind(all_baseline_files(idx,1).name, pattern);
    if notTrueMAT
        ; % do nothing for now.  maybe delete these files in future. 
    else
        load(all_baseline_files(idx,1).name);
        F0s = [F_0; F0s]; 
    end
end
clear idx; 
clear pattern;

% Assert that all trials being plotted together were acquired with the same
% time points (same number of frames and length of period).
if size(tPts,1) > 1
    assert(max(max(tPts) - min(tPts)) == 0,...
        'Expect all timePts to be the same') ;
end

% Calculate area under curve for each trial, starting at time 3s.
% Get index of odor onset frame: 
for idx = 1:length(tPts)
    if ~(tPts(1,idx) <onTime)
        onIdx = idx;
        break;
    end
end
% Get index of odor offset frame: 
for idx = onIdx:length(tPts)
    if (tPts(1,idx) > offTime)
        offIdx = idx;
        break;
    end
end

% postOnset = []; % init
AOCsum = nan(size(all_baseline_files));
Atrapz = nan(size(all_baseline_files));
dfPeak = nan(size(all_baseline_files)); % peak of each trial
for trial =1:size(tPts,1)
    if trial == 1
        postOnset = nan(size(dFoF,1),size(dFoF,2)+1-onIdx);
    end
    temp = dFoF(trial,onIdx:end);
    odorONdFoF = dFoF(trial,onIdx:offIdx);
%     postOnset = [temp; postOnset];
    postOnset(trial,:) = temp;
    AOCsum(trial) = sum(postOnset(trial,:)); 
    Atrapz(trial) = trapz(postOnset(trial,:));
    AOdortrapz(trial) = trapz(odorONdFoF); % between timeOn and timeOff
    dfPeak(trial) = max(dFoF(trial,:)); % peak in entire trial    
end
%% Plot GCaMP (for each trial and on average) versus time. 
% framePeriod = 0.818748; % Frame period : seconds/frame. 
% framePeriod fetched automatically from config file, stored as framePerVal. 
% timePoints = [ 1:maxTime ] * framePeriod; 
fig = figure;
xvalsOn = [3,3];
xvalsOff = [5,5];
% plot(timePoints, CaAvg);
% plot(timePoints, delFoverF);
subplot(1,2,1); % trial-by-trial plotting. 
hold all
evokedTrialSums = zeros(size(dFoF, 1),1); % initialize 
for idx = 1:max(size(dFoF, 1))
    plot(tPts(1,:), dFoF(idx,:),'LineWidth',3);
    %
    ind = 0; 
    for i = 1:length(timePoints)
        if (timePoints(i) > onTime) && (timePoints(i) < offTime)
            ind = ind +1;
            evokedTrialSums(idx) = evokedTrialSums(idx) + dFoF(idx, i); 
            denom = ind; 
        end
    end
end
% set limits of x axis
set(gca,'xlim', x_limits);

evokedTrialAvgs = evokedTrialSums/denom; 
%
%legend('show')
legend(stringsList) % "creates a legend in the current axes and uses the entries in strings to label each set of data. Specify strings as either a cell array of strings or a matrix of strings."
% mean(A,1) is a row vector containing the mean value of each column.
avg_dFoF = mean(dFoF,1);
max_dFoF = max(dFoF);
min_dFoF = min(dFoF);

% % Calculate and set max value for y axis. 
% % dFoF_max = max(max(dFoF));
% % avg_dFoF_max = max(max(avg_dFoF));
% % set_y_max = max(dFoF_max, avg_dFoF_max);

title('deltaF/F vs Time')
% Annotate the point (3, delFoverF(3))
% text(3, delFoverF(3),'\leftarrow (delFoverF(3sec))',...
%      'HorizontalAlignment','left')
% title('GCaMP Average vs Time')

% % Add a line at timepoint = 3 sec
    yvals = get(gca, 'ylim');
    max_y_vals = yvals
    
    lineOn = plot(xvalsOn, yvals);
    lineOff = plot(xvalsOff, yvals);
    % Make the new vertical line at 3s a dash-dotted black line. 
    set(lineOn,'Color','k','LineWidth', 1, 'LineStyle','-.')
    set(lineOff,'Color','k','LineWidth', 1, 'LineStyle','-.')
hold off
hLegend1 = legend({'Trial1','Trial2','Trial3','Trial4','Trial5'}); 
hYLabel = ylabel('dF/F_0');

% Also plot the mean of each column (timePt):
subplot(1,2,2);
avg_dFoF_plot = plot(tPts, avg_dFoF,'Color','b','LineWidth',3);

hold on
hTitle2 = title('Average deltaF/F vs Time');

% Annotate the point (3, delFoverF(3))
% % % text(3, delFoverF(3),'\leftarrow (delFoverF(3sec))',...
% % %      'HorizontalAlignment','left');

% % % Add max and min curves.
% % maxCurve = plot(tPts(1,:), max_dFoF);
% % set(maxCurve, 'Color', 'blue', 'LineWidth', 1, 'LineStyle', '--');
% % minCurve = plot(tPts(1,:), min_dFoF);
% % set(minCurve, 'Color', 'blue', 'LineWidth', 1, 'LineStyle', '--');
%
% How to fill in space between curves (not working properly):
% curveDiffVec = [min_dFoF, max_dFoF];
% xCurve = [timePoints,timePoints];
% fill(xCurve, curveDiffVec, 'b'); % In the style of: http://www.mathworks.com/matlabcentral/answers/14005
% Add a line at timepoint = 3 sec
    set(gca, 'ylim', max_y_vals);
    yvals = get(gca, 'ylim');
%     xvalsOn = [3,3];
    lineOn = plot(xvalsOn, yvals);
%     xvalsOff = [4,4];
    lineOff = plot(xvalsOff, yvals);
    % Make the new vertical line at 3s a dash-dotted black line. 
    set(lineOn,'Color','k','LineWidth', 1, 'LineStyle','-.')
     % Make the new vertical line at 4s a dash-dotted black line.
    set(lineOff,'Color','k','LineWidth', 1, 'LineStyle','-.')
    %legend
hold off
set(gca,'xlim', x_limits);
hLegend2 = legend('Average');
hXLabel2 = xlabel('Time (Seconds)'); 
hYLabel = ylabel('dF/F_0'); 

%% Single-pane figure with trials and average trace. 
fig = figure;
xvalsOn = [3,3];
xvalsOff = [5,5]; 
hold on
rectColor = 0.85*[1 1 1]; 
rectangle('Position',[onTime,-0.08,offTime-onTime,0.05],'FaceColor',rectColor,'EdgeColor',rectColor,...
    'LineWidth',3)
evokedTrialSums = zeros(size(dFoF, 1),1); % initialize 
for idx = 1:max(size(dFoF, 1))
    plot(tPts(1,:), dFoF(idx,:),'LineWidth',3);
    ind = 0; 
    for i = 1:length(timePoints)
        if (timePoints(i) > onTime) && (timePoints(i) < offTime)
            ind = ind +1;
            evokedTrialSums(idx) = evokedTrialSums(idx) + dFoF(idx, i); 
            denom = ind; 
        end
    end
end
plot(tPts, avg_dFoF,'Color','k','LineWidth',3);
set(gca,'xlim', x_limits); % set limits of x axis
evokedTrialAvgs = evokedTrialSums/denom; 
legend(stringsList); 
avg_dFoF = mean(dFoF,1);
max_dFoF = max(dFoF);
min_dFoF = min(dFoF);
sem_dFoF = nan(size(max_dFoF)); 
for timei = 1:length(dFoF)
    sem_dFoF(1,timei) = get_sem(dFoF(:,timei)); 
end

hold off
hTitle = title('Odor Response Example'); 
hXLabel = xlabel('Time (Seconds)'); 
hYLabel = ylabel('dF/F_0');
hText = text(3.08,-.05,'Odor On','Color','k','FontSize',12);
hLegend = legend({'Trial1','Trial2','Trial3','Trial4','Trial5','Average'});
% Make it prettier
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel, hText], ...
    'FontName'   , 'Helvetica');
set([hLegend, gca]             , ...
    'FontSize'   , 12           );
set([hXLabel, hYLabel, hText]  , ...
    'FontSize'   , 16          );
set( hTitle                    , ...
    'FontSize'   , 16          , ...
    'FontWeight' , 'bold'      );
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'      , 0:0.6:0.6, ...
  'LineWidth'   , 1         );

try
    print -depsc '/Users/katherineshakman/Dropbox/Axel Lab Dropbox/JClub and LM Papers and Notes/MyPresentations/Figures/ExampleOdor.eps'
catch
    print -depsc '/Users/katieshak/Dropbox/Axel Lab Dropbox/JClub and LM Papers and Notes/MyPresentations/Figures/ExampleOdor.eps'
end
%% Plot with mean and sem shaded error bars
figure; 
hold on
rectColor = 0.85*[1 1 1]; 
rectangle('Position',[onTime,-0.08,offTime-onTime,0.05],'FaceColor',rectColor,'EdgeColor',rectColor,...
    'LineWidth',3)

plot(tPts, avg_dFoF,'Color','k','LineWidth',3); 
shadedErrorBar(tPts(1,:), avg_dFoF, sem_dFoF); 

set(gca,'xlim', x_limits);
set(gca,'ylim', max_y_vals); 
hLegend2 = legend('Average');
hXLabel2 = xlabel('Time (Seconds)'); 
hYLabel = ylabel('dF/F_0'); 

hTitle = title('Odor Response Example'); 
hXLabel = xlabel('Time (Seconds)'); 
hYLabel = ylabel('dF/F_0');
hText = text(3.08,-.05,'Odor On','Color','k','FontSize',12);
hLegend = legend({'Average +/- SEM'});
% Make it prettier
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel, hText], ...
    'FontName'   , 'Helvetica');
set([hLegend, gca]             , ...
    'FontSize'   , 12           );
set([hXLabel, hYLabel, hText]  , ...
    'FontSize'   , 16          );
set( hTitle                    , ...
    'FontSize'   , 16          , ...
    'FontWeight' , 'bold'      );
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'      , 0:0.6:0.6, ...
  'LineWidth'   , 1         );

try
    print -depsc '/Users/katherineshakman/Dropbox/Axel Lab Dropbox/JClub and LM Papers and Notes/MyPresentations/Figures/ExampleOdorSEM.eps'
catch
    print -depsc '/Users/katieshak/Dropbox/Axel Lab Dropbox/JClub and LM Papers and Notes/MyPresentations/Figures/ExampleOdorSEM.eps'
end
%% Save figure files and evokedAvg value
dateStart = strfind(imDir, '201'); 
dateEnd = dateStart-1 + min(strfind(imDir(dateStart:end), '/'))-1; 
date = imDir(dateStart:dateEnd);
trialStart = max(strfind(imDir, '/'))+1; 
trialID = imDir(trialStart:end); 
figName = [date,'_', trialID];
%savePath = currentDirectory(...
            %strfind(currentDirectory, 'Raw_Un')+16:strfind(...
            %currentDirectory, 'TSeries')-1);
%figName = savePath; 
saveas(fig, figName, 'fig'); 
saveas(fig, figName, 'png'); 

ind = 0; 

for i = 1:length(timePoints)
    if (timePoints(i) > 3) && (timePoints(i) < 5)
        ind = ind +1;
        evokedTimes(ind) = timePoints(i);
        evokedResp(ind) = avg_dFoF(i);
    end
end

% convert currentDirectory to currentDirName
currentDirName = strrep(currentDirectory, '/', '_')  % Print to console.

evokedAvg = mean(evokedResp); 
evokedAvg_savefile = [figName, '_evokedAvg']; 
save(evokedAvg_savefile, 'evokedAvg'); 
load(evokedAvg_savefile); 
%
evokedPeak = max(evokedResp); 
evokedPeak_savefile = [figName, '_evokedPeak']; 
save(evokedPeak_savefile, 'evokedPeak'); 
load(evokedPeak_savefile); 
%
avgBaseline = mean(F0s); 
avgBaseline_savefile = [figName, '_avgBaseline'];
save(avgBaseline_savefile, 'avgBaseline'); 
load(avgBaseline_savefile); 
% 
tPts_savefile = [figName, '_tPts'];
save(tPts_savefile,'tPts'); % save timePoints 
%

% %display('evokedAvg ='); 
% evokedAvg
% %display('evokedPeak ='); 
% evokedPeak
% %
% avgBaseline
% %
% evokedTrialAvgs % display values

avg_dFoF_savefile = [figName, '_avg_dFoF'];
save(avg_dFoF_savefile, 'avg_dFoF'); 

Atrapz_savefile = [figName, 'AOC_trapz'];
cd('..'); % up one level
save(Atrapz_savefile, 'Atrapz');

OdorAtrapz_savefile = [figName, 'AOdortrapz'];
save(OdorAtrapz_savefile, 'AOdortrapz');

avgAOdortrapz = mean(AOdortrapz); 
avg_OdorAtrapz_savefile = [figName, 'avg_AOCOdor_trapz'];
save(avg_OdorAtrapz_savefile, 'avgAOdortrapz');

avgAtrapz = mean(Atrapz); 
avg_Atrapz_savefile = [figName, 'avg_AOC_trapz'];
save(avg_Atrapz_savefile, 'avgAtrapz');

dfPeak_savefile = [figName, '_trial_dfPeak'];
% cd('..'); % up one level
save(dfPeak_savefile, 'dfPeak');

avgdfPeak = mean(dfPeak);
avg_dfPeak_savefile = [figName, '_avg_dfPeak'];
save(avg_dfPeak_savefile, 'avgdfPeak');

% % % save a copy in the GroupedData folder)
% % savepath = '/Volumes/KATIELAB1/Data_Analysis/GroupData/' ;
% % cd(savepath);
% % save(evokedAvg_savefile, 'evokedAvg');
% % save(avg_dFoF_savefile, 'avg_dFoF', 'tPts'); 

% save(fullfile(savepath,evokedAvg_savefile,'evokedAvg'));
% save(fullfile(savepath,avg_dFoF_savefile,'avg_dFoF'));
%save(avg_dFoF_savefile, 'avg_dFoF'); 

cd(StartDir); 