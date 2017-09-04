%% Plot Only Fluor Analysis with Z-score
% Plot from previously stored deltaF/F data and timepoints.
clear; close all;

%function [] = PlotAll_DelFoverF_ofType_zScore()

StartDir = pwd;

% Variables user can change:
x_limits = [0,10];
y_line_limits = [-1,1];
START_BASELINE = 2;
END_BASELINE = 3;
END_ODOR = 5;

%% Load data files from several user-selected directories.

imDir = pwd;
currentDirectory = pwd; % For use in save paths.

all_delFoverF_files = dir('*delFoverF.mat');
all_delFoverF_files = flipud(all_delFoverF_files);
all_timePoints_files = dir('*timePoints.mat');
all_timePoints_files = flipud(all_timePoints_files);
all_baseline_files = dir('*baselineF0.mat');
all_baseline_files = flipud(all_baseline_files);

dFoF = [];
tPts = {}; % cell
tPtsCount = 0; % initialize
longestTimePts = 0:0.2:10.2; % initialize
F0s = [];
stringsList = []; % initialize
dFoF = {};
dFcount = 0; % initialize
% dFoF = zeros(1, max(size(all_delFoverF_files))); % Initialize
% tPts = zeros(1, max(size(all_timePoints_files))); % Initialize
pattern = '._' ;

for trial = 1:max(size(all_timePoints_files))
    notTrueMAT = strfind(all_timePoints_files(trial,1).name, pattern);
    if notTrueMAT
        ; % do nothing for now.  maybe delete these files in future.
    else
        load(all_timePoints_files(trial,1).name);
        tPtsCount = tPtsCount + 1; % iterate
        tPts{tPtsCount,1} = timePoints;
        if length(timePoints) > length(longestTimePts)
            longestTimePts = timePoints;
        end
    end
end

for trial = 1:max(size(all_delFoverF_files))
    %     new_delFoverF = load(all_delFoverF_files(idx,1).name);
    try
        notTrueMAT = strfind(all_delFoverF_files(trial,1).name, pattern);
    catch
        notTrueMAT = NaN;
        notTrueMAT = false;
    end
    if notTrueMAT
        ; % do nothing for now.  maybe delete these files in future.
    else
        % Load dF/F data
        load(all_delFoverF_files(trial,1).name);
        dFcount = dFcount + 1; % iterate
        dFoF{dFcount,1} = delFoverF;
        stringsList = [all_delFoverF_files(trial,1).name(9:25); stringsList];
        % used only part of the name that contains date and tseries number
    end
end

% Interpolate dF/F data
% Want to interpolate so we can look at all of these traces
% together.
disp('Interpolating...');
preInterp_dFoF = dFoF;
dFoF = cell(size(tPts));
for trialIdx = 1:length(tPts)
    tempInterp = interp1(tPts{trialIdx},preInterp_dFoF{trialIdx},longestTimePts);
    dFoF{trialIdx} = tempInterp;
end


stringsList = cellstr(stringsList); % make cell array list of strings

for trial = 1:max(size(all_baseline_files)) % also get baselines
    notTrueMAT = strfind(all_baseline_files(trial,1).name, pattern);
    if notTrueMAT
        ; % do nothing for now.  maybe delete these files in future.
    else
        load(all_baseline_files(trial,1).name);
        F0s = [F_0; F0s];
    end
end

clear pattern;

%% Check data and Calc z-scores

% Calculate area under curve for each trial, starting at time 3s.
for trial = 1:length(timePoints)
    if ~(timePoints(1,trial) <END_BASELINE)
        onIdx = trial;
        break;
    end
end
% postOnset = []; % init
AOCsum = nan(size(all_baseline_files));
Atrapz = nan(size(all_baseline_files));
dfPeak = nan(size(all_baseline_files)); % peak of each trial
% raw = []; % init
for trial =1:length(tPts)
    assert(size(tPts,2) < size(tPts,1)); % if fails, probably have huge num trials or too few frames
    if trial == 1
        postOnset = nan(size(dFoF{1},1),size(dFoF{1},2)+1-onIdx);
    end
    temp = dFoF{trial}(onIdx:end);
    %     postOnset = [temp; postOnset];
    postOnset(trial,:) = temp;
    AOCsum(trial) = sum(postOnset(trial,:));
    Atrapz(trial) = trapz(postOnset(trial,:));
    dfPeak(trial) = max(dFoF{trial}); % peak in entire trial
    % Get back the raw traces:
    raw(trial,:) = (dFoF{trial}+F0s(trial))*F0s(trial);
    % Get mean and sigma during baseline (pre-stimulus) for the trace:
    preOdorInd = longestTimePts<END_BASELINE; % gets recalculated each trial
    preOdorRaw = raw(trial,: ) .* preOdorInd; % gets recalculated each trial
    postStartBaseInd = longestTimePts>START_BASELINE;
    baselnInd = preOdorInd.*postStartBaseInd;
    baselnRaw = raw( trial, : ) .* baselnInd;
    raw_preOdor_mean(trial) = nanmean(preOdorRaw); % single value each trial
    raw_preOdor_sigma(trial) = nanstd(preOdorRaw); % single value each trial
    raw_baseln_mean(trial) = nanmean(baselnRaw);
    raw_baseln_sigma(trial) = nanstd(baselnRaw);
    % Get zscores for each trace:
    zscoresLongBaseln(trial,:) = ( raw(trial,:)-raw_preOdor_mean(trial) )/raw_preOdor_sigma(trial) ;
    zscores(trial,:) = ( raw(trial,:)-raw_baseln_mean(trial) ) / raw_baseln_sigma(trial) ;
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
evokedTrialSums = zeros(length(dFoF),1); % initialize
for trial = 1:length(dFoF)
    plot(longestTimePts, dFoF{trial});
    %
    ind = 0;
    for i = 1:length(timePoints)
        if (timePoints(i) > END_BASELINE) && (timePoints(i) < END_ODOR)
            ind = ind +1;
            evokedTrialSums(trial) = evokedTrialSums(trial) + dFoF{trial}(i);
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
dFoFmat = cell2mat(dFoF);
avg_dFoF = mean(dFoFmat,1);
max_dFoF = max(dFoFmat);
min_dFoF = min(dFoFmat);

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

try
    lineOn = plot(xvalsOn, y_line_limits);
    lineOff = plot(xvalsOff, y_line_limits);
catch
    lineOn = plot(xvalsOn, yvals);
    lineOff = plot(xvalsOff, yvals);
end
% Make the new vertical line at 3s a dash-dotted black line.
set(lineOn,'Color','k','LineWidth', 1, 'LineStyle','-.')
set(lineOff,'Color','k','LineWidth', 1, 'LineStyle','-.')
hold off
ylim(max_y_vals); % rescale the y-axis in case the vertical lines are too tall relative to traces

% Also plot the mean of each column (timePt):

subplot(1,2,2);
avg_dFoF_plot = plot(longestTimePts, avg_dFoF);
set(avg_dFoF_plot, 'Color', 'blue', 'LineWidth', 1)

hold on
title('Average deltaF/F vs Time')

% Annotate the point (3, delFoverF(3))
% % % text(3, delFoverF(3),'\leftarrow (delFoverF(3sec))',...
% % %      'HorizontalAlignment','left');

% Add a line at timepoint = 3 sec
set(gca, 'ylim', max_y_vals);
yvals = get(gca, 'ylim');
lineOn = plot(xvalsOn, yvals);
lineOff = plot(xvalsOff, yvals);
% Make the new vertical line at 3s a dash-dotted black line.
set(lineOn,'Color','k','LineWidth', 1, 'LineStyle','-.')
% Make the new vertical line at 4s a dash-dotted black line.
set(lineOff,'Color','k','LineWidth', 1, 'LineStyle','-.')
%legend
hold off
set(gca,'xlim', x_limits);

%% Plot z-scores(for each trial) versus time.
zfig = figure;
xvalsOn = [3,3];
xvalsOff = [5,5];
% subplot(1,2,1); % trial-by-trial plotting.
hold all
evokedTrialZSums = zeros(size(zscores, 1),1); % initialize
for trial = 1:length(dFoF)
    plot(longestTimePts, zscores(trial,:));
    %
    ind = 0;
    for i = 1:length(timePoints)
        if (timePoints(i) > END_BASELINE) && (timePoints(i) < END_ODOR)
            ind = ind +1; % counting frames during stimulus
            evokedTrialZSums(trial) = evokedTrialZSums(trial) + zscores(trial, i);
            denom = ind;
        end
    end
end
set(gca,'xlim', x_limits); % set limits of x axis

evokedTrialZAvgs = evokedTrialZSums/denom;

legend(stringsList) % "creates a legend in the current axes and uses the entries in strings to label each set of data. Specify strings as either a cell array of strings or a matrix of strings."
% mean(A,1) is a row vector containing the mean value of each column.
avg_dFoF = mean(cell2mat(dFoF),1);
max_dFoF = max(cell2mat(dFoF));
min_dFoF = min(cell2mat(dFoF));

title('z-score vs Time (mean and sd calculated during baseline period)')
% Annotate the point (3, delFoverF(3))
% text(3, delFoverF(3),'\leftarrow (delFoverF(3sec))',...
%      'HorizontalAlignment','left')
% title('GCaMP Average vs Time')

% % Add a line at timepoint = 3 sec
yvals = get(gca, 'ylim');
max_y_vals = yvals

try
    lineOn = plot(xvalsOn, max_y_vals);
    lineOff = plot(xvalsOff, max_y_vals);
catch
    lineOn = plot(xvalsOn, yvals);
    lineOff = plot(xvalsOff, yvals);
end
% Make the new vertical line at 3s a dash-dotted black line.
set(lineOn,'Color','k','LineWidth', 1, 'LineStyle','-.')
set(lineOff,'Color','k','LineWidth', 1, 'LineStyle','-.')
hold off
ylim(max_y_vals); % rescale the y-axis in case the vertical lines are too tall relative to traces

%% Plot full length baseline z-scores avg versus time.
zLongFig = figure;
xvalsOn = [END_BASELINE,END_BASELINE];
xvalsOff = [END_ODOR,END_ODOR];
zscoresAvgLong = mean(zscoresLongBaseln,1);

hold all
plot(longestTimePts, zscoresAvgLong);
set(gca,'xlim', x_limits); % set limits of x axis
legend('Average') % "creates a legend in the current axes and uses the entries in strings to label each set of data. Specify strings as either a cell array of strings or a matrix of strings."
title('Long Baseline z-score Avg vs Time (mean and sd calculated during full baseline period)')

% mean(A,1) is a row vector containing the mean value of each column.
avg_dFoF = mean(cell2mat(dFoF),1);
max_dFoF = max(cell2mat(dFoF));
min_dFoF = min(cell2mat(dFoF));

% % Add a line at timepoint = 3 sec
yvals = get(gca, 'ylim');
max_y_vals = yvals;

try
    lineOn = plot(xvalsOn, max_y_vals);
    lineOff = plot(xvalsOff, max_y_vals);
catch
    lineOn = plot(xvalsOn, yvals);
    lineOff = plot(xvalsOff, yvals);
end
% Make the new vertical line at 3s a dash-dotted black line.
set(lineOn,'Color','k','LineWidth', 1, 'LineStyle','-.')
set(lineOff,'Color','k','LineWidth', 1, 'LineStyle','-.')
hold off
ylim(max_y_vals); % rescale the y-axis in case the vertical lines are too tall relative to traces

%% Make a single figure with the mean and SEM error bars for all the traces in this group/folder:
avgdfEachTimepoint = mean(cell2mat(dFoF),1);
ntrials = length(dFoF)
SDdfEachTimepoint = std(cell2mat(dFoF),1);
SEMdfEachTimepoint = SDdfEachTimepoint./sqrt(ntrials);
figure;
plot(longestTimePts,avgdfEachTimepoint);
hold on
errorbar(longestTimePts,avgdfEachTimepoint,SEMdfEachTimepoint);
hold off
%% Make a single-panel figure with the lines and the average superimposed in black line:
% newFig = figure;
% new_evokedTrialSums = zeros(size(dFoF, 1),1); % initialize
% hold on
% for idx = 1:max(size(dFoF, 1))
%     plot(tPts(1,:), dFoF(idx,:));
%     ind = 0;
%     for i = 1:length(timePoints)
%         if (timePoints(i) > 3) && (timePoints(i) < 5)
%             ind = ind +1;
%             new_evokedTrialSums(idx) = new_evokedTrialSums(idx) + dFoF(idx, i);
%             denom = ind;
%         end
%     end
% end
% % set limits of x axis
% set(gca,'xlim', x_limits);
%     yvals = get(gca, 'ylim');
%     max_y_vals = yvals
%     lineOn = plot(xvalsOn, yvals);
%     lineOff = plot(xvalsOff, yvals);
%     set(lineOn,'Color','k','LineWidth', 1, 'LineStyle','-.')    % Make the new vertical line at 3s a dash-dotted black line.
%     set(lineOff,'Color','k','LineWidth', 1, 'LineStyle','-.')
% ylim(max_y_vals); % rescale the y-axis in case the vertical lines are too tall relative to traces
% avg_dFoF_plot = plot(tPts, avg_dFoF);
% hold off
% xlabel('Seconds'); ylabel('dF/F');
%
% a=findobj(gcf); % get the handles associated with the current figure
%
% allaxes=findall(a,'Type','axes');
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
%
% set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,...
%     'FontSize',14);
% set(alllines,'Linewidth',4);
% set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14);
% % info on setting all axes, lines, and text properties from: http://matlab.cheme.cmu.edu/2011/08/01/plot-customizations-modifying-line-text-and-figure-properties/
% set(avg_dFoF_plot, 'Color', 'black', 'LineWidth', 5)

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