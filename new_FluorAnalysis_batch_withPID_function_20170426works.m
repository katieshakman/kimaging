function newFluorAnalysis_batch_withPID_func(ROI, roiName, inputFolder)

% %% For testing: 
% inputFolder = pwd; 
% %'/Users/katieshak/Desktop/Data_Analysis_Temp/Raw_Unprocessed/2015_10_08_testing/br1_58B_6fA/p1_hemi1';
% roiName = 'roiName'; 

%%% Fluorescence Analysis
% Code to process GCaMP calcium imaging Tseries (from 2-photon).
% What it does/will do:
% Reads image, allows selection of an ROI (previously designed), gets fluorescence information
% for the ROI over time, saves a chart of the deltaF/F values, and makes a
% plot of deltaF/F for the ROI over time.

% Imaging channels:
% Ch1 : red PMT 1 HV (RFP)
% Ch2 : green gaasp (GCaMP or GFP)
% Ch3 : dodt (transmitted light)

% clear all; close all;

%%% Option to add ROIs interactively. 
%useROI = 0; % set to 0 for no ROI, or to 1 for one or more ROIs.  
%%% Option to gather PID data. 

%%% Option to add ROIs interactively.
%useROI = 1; % set to 0 for no ROI, or to 1 for one or more ROIs.
%%% Option to gather PID data.
%usePID = 0; % set to 0 for ignoring PID data, or to 1 to get/plot PID data.
usePID = 0; 

% Figure out from input if using ROI already: 
if exist('ROI','var') == 0 % if no ROI was specified
    useROI = 0; % then set useROI to 0
elseif isempty(ROI) % if there is an empty ROI variable
    useROI = 0; % then set useROI to 0
else % if there is a non-empty ROI variable
    useROI = 1; % then set useROI to 1
end

%% Set up the folders
myFolder = pwd;
stimFolder = strrep(inputFolder,'Raw_Unprocessed','Data_Analysis'); % modifiedStr = strrep(origStr, oldSubstr, newSubstr)
disp('Analysis folder is: '); disp(stimFolder); % the analysis folder should be a short stimulus name, like 2-Oct
exptFolder = strrep(myFolder,inputFolder,'');

%% Compile list of Tseries directories in current folder.
TS_listing = dir('Tseries*');
if isempty(TS_listing)
    TS_listing = dir('TSeries*'); % Try again (this sometimes works). 
    if isempty(TS_listing)
        display('No TSeries found in directory.')
    end
end

%%% Read the image (and get an ROI, if desired).
for TS_idx = 1:length(TS_listing)
    %close all;
    % clear all variables except myFolder, TS_listing, and TS_idx:
    % clearvars -except myFolder TS_listing TS_idx LongSavePath useROI ROI roiName usePID inputFolder stimFolder exptFolder;
    %     clear all; close all;
    %     imDir = uigetdir(); % Allow user to specify the directory to read images from.
    mySubfolder = TS_listing(TS_idx).name;
    imDir = fullfile(myFolder, mySubfolder);
    cd(imDir);
    tifIm = dir('*.tif');
    %im1 = imread(tifIm(1).name); % Below: changed this to the average image for the stack.
    cd ../; % Move to containing folder (up one level).
    currFolderName = pwd; % Get name of containing folder.
    cd(imDir);
    
    getAllChAvg = 1;
    if getAllChAvg == 1
        % Get the average image for the stack in all channels:
        imSize = size(imread(tifIm(1).name));
        runningTotIm = uint16(zeros(imSize)); % Initialize
        for idx = 1:length(tifIm)
            newTifIm = imread(tifIm(idx).name);
            lastIdx = idx; % Save the last idx value for the average.
            runningTotIm = runningTotIm + newTifIm;
        end
        % avgIm = zeros(512,512); % Initialize the average image variable.
        avgIm = runningTotIm./length(tifIm);
    end
    % Make ROI if desired
    if useROI == 1
        if ~exist('ROI') || isempty('ROI')
            disp('Please draw an ROI on the figure and double-click when done.')
            ROI = roipoly(runningTotIm);
            ROI = uint16(ROI);
            roiName = 'new';
        elseif isnan(ROI);
            ROI = uint16(ones(size(avgIm))); % Use entire average image size for ROI if NaN ROI was passed.
            roiName = 'full';
        else
            figure;
            subplot(1,2,1); imshow(runningTotIm);
            subplot(1,2,2); imshow(ROI.*(runningTotIm));
            % TODO: Add a dialog to confirm or adjust ROI.
        end
        % Mask average image with ROI.
        avgIm = ROI.*avgIm;
    elseif useROI == 2
        % load ROI
    else
        ROI = ones(size(avgIm));
    end
    
    % While in imDir, read and save the PID info.
    if usePID == 1
        PID_list = dir('*.csv');
        [PID_fig, PID_figName] = PID_Plot_func(PID_list(1,1).name);
    end
    
    %% Get the value of framePer for this image stack. Modified From RosettaCode
    % "This is defined as a function, parameters are returned as part
    % of a struct. When the first line, and the assignment to return
    % values are removed, it is a script that stores the parameters in
    % the local workspace.
    % Source:
    % http://rosettacode.org/wiki/Read_a_configuration_file#MATLAB_.2F_Octave
    
    % configfile = dir('*.cfg'); % Changed from uigetfile to dir
    %c = dir('*.cfg');
    % %     c = dir('*.xml');
    % %     configfile = c.name; %c(1,1).name;
    % %     fid = fopen(configfile);
    % %     if fid<0, error('cannot open file %s\n',a); end;
    % %
    % %     while ~feof(fid)
    % %         line = strtrim(fgetl(fid));
    % %         framePerLocation = strfind(line, 'framePeriod=');
    % %         if ~isempty(framePerLocation)
    % %             framePerLine = line;
    % %             %     else isempty(line) || all(isspace(line)) || strncmp(line,'#',1) || strncmp(line,';',1),
    % %             %         ; % no operation
    % %         else
    % %             framePerLine = framePerLocation;
    % %         end;
    % %     end;
    % %     fclose(fid);
    % %     % whos,     % shows the parameter in the local workspace
    
    c = dir('TS*.xml');
    myXML = c.name;
    mytext = fileread(myXML);
    framePerLoc = strfind(mytext, 'framePeriod" value="');
    framePerStrLen = length('framePeriod" value="');
    framePerEnd = framePerLoc + framePerStrLen + 10;  % can be adjusted to auto-length later
    framePerVal = str2double(mytext(framePerLoc+framePerStrLen:framePerEnd));
    % Now copy the repetitionPeriod value out of the line we just found:
    %repPerLoc = strfind(framePerLine, 'repetitionPeriod="');
    % repetitionPeriod seems to be framePeriod (probably not scanLinePeriod) in new prairieView
    % repPerStrLength = length('repetitionPeriod="'); % 18 for 'repetitionPeriod="'
    %repPerStrLength = length('framePeriod" value="');
    %endLoc = strfind(framePerLine, 'triggerSyncExperiment=') - 3; % changed: from ? ? 2? to ? ? 3?.
    %framePerVal = str2double(mytext(framePerLoc:framePerEnd)); % Maybe change this so all digits are preserved, instead of truncated?
    % framePerVal = str2ndouble(framePerLine(valueLoc+7:end-5));
    
    getCh1Avg = 0; % Set this value to 1 to get the Ch1 Avg Img.
    if getCh1Avg == 1
        % Get the average image for the stack in only Ch1:
        frameSize = imSize; % Auto-adjusts for size of frames in the Tseries.
        runningTotIm = uint16(zeros(frameSize)); % Initialize
        pattern = 'Ch1'; % Pattern to search for in image name string.
        for idx = 1:length(tifIm)
            % Locate "Channel 2" in current image name:
            k = strfind(tifIm(idx).name, pattern);  % Returns empty var if pattern not found.
            if ~isempty(k) % Read im if k is NOT empty (Read Channel 2 images).
                newTifIm = imread(tifIm(idx).name);
            end
            runningTotIm = runningTotIm + newTifIm;
        end
        avgIm = runningTotIm./(length(tifIm)/3); % TODO: Find the actual # of Ch2 Ims.
        clear runningTotIm;
        clear k;
    end
    
    % Make the figure(s).
    % reply = input('Do you wish to select an ROI? y/n (Default is n)', 's');
    reply = 'n'; % User may select an ROI when reply = y;
    % figure;
    divImBy = 2; % Daisuke initially suggested 16 for this value.
    
    % subplot(1,2,1); imshow(uint8(avgIm))
    % subplot(1,2,2); imshow(uint8(avgIm)./divImBy)
    % title('figure1')
    % if isempty(reply)
    %     reply = n;
    % elseif strcmp(reply,'y')
    %     x = impoly;
    %     pos = getPosition(x);
    %     mask = createMask(x);
    %     figure; imshow(mask)
    %     title('figure2')
    % %     roiv1 = double(im1).*mask;
    %     roiv1 = double(avgIm).*mask;
    %     figure; imshow(uint8(roiv1))
    % %     figure; imshow(uint8(roiv1./16))
    %     title('figure3')
    % end
    
    %%% Get GCaMP average at each time point; convert to delF/F.
    % TODO: Add in ROI restriction.
    % Do this only for the part of tifIm that is Ch1 (GCaMP)
    
    %Ch1tifIm = %all of tifIm for which tifIm(1).name contains Ch1
    % TODO: Check that there are Ch1, Ch2, and Ch3 images
    % before dividing by 3 (if only Ch1 and Ch2, divide by 2).
    % % % maxTime = length(tifIm)/3; % Last timepoint/slice,
    % div by 3 since only first 1/3 of the Tseries files will be Ch1.
    % TODO: Change to Ch1tifIm instead of tifIm, so we are sure we only use Ch1
    % images.
    
    tifIm_names = cell(length(tifIm),1);
    for idx = 1:length(tifIm)
        tifIm_names{idx} = tifIm(idx).name;
    end
    Ch1_indices = strfind(tifIm_names, 'Ch1');
    tifImCh1_indices = cell2mat(Ch1_indices);
    % tifImCh1_indices = tifImCh1_indices * ones(length(tifImCh1_indices));
    maxTime = length(tifImCh1_indices);
    CaAvgDenom = 512*512; % Total number of pixels per slice.
    % TODO: Change the above (CaAvgDenom) to get the true image (or ROI) size.
    
    for timeIdx = 1:maxTime
        currTif = imread(tifIm(timeIdx).name) ;
        % Mask all images in stack and avg image with ROI (if using):
        if useROI
            currTif = ROI.*currTif;
        end
        CaAvg(timeIdx) = sum(sum( currTif )) / CaAvgDenom; % Avg Signal in frame normalized to the image size.
    end
    timePt3sec = 3/framePerVal; % 3sec*(1/(sec/frame)) = 3sec/framePerVal
    if timePt3sec < 1 || isnan(framePerVal)
        F_0 = mean(CaAvg(1:1)); % Average fluorescence in first frame.
    else
        F_0 = mean(CaAvg(1:floor(timePt3sec))); % Average fluorescence up to 3 sec (takes floor of timePt3sec).
    end
    
    delF = CaAvg - F_0*ones(size(CaAvg));
    delFoverF = delF/F_0;
    timePoints = [ 1:maxTime ] .* framePerVal;
    % Above: time is in sec, framePerVal is in sec/frame --> divide by framePer
    % allTimePoints = [framePerVal:maxTime];
    
    %% Make a heatmap from 3 sec to 5 sec
    % get image after stimulus has been on for 1 second (this should be
    % time = 4sec)
    timePt5sec = 5/framePerVal; % to get frame # at which time is 4 sec
    sec5TifIm = imread(tifIm(ceil(timePt5sec)).name);
    % Get average baseline image:
    % get baseline image, then divide by that image (so each point in image is converted
    % to its own delF/F0 value).
    baseImSum = im2double(imread(tifIm(1).name)); % Inialize to first image.
    for imIdx=2:ceil(timePt3sec) % Read in baseline images.
        tempIm = imread(tifIm(imIdx).name);
        baseImSum = baseImSum + im2double(tempIm);
    end
    meanBaseIm = baseImSum/ceil(timePt3sec);
    clear tempIm;
    % Get average stimulus image:
    stimImSum = im2double(imread(tifIm(ceil(timePt3sec)).name)); % Inialize to first stim image.
    for imIdx = 2:ceil(timePt5sec); % Read in stimulus images.
        tempIm = imread(tifIm(imIdx).name);
        stimImSum = stimImSum + im2double(tempIm);
    end
    meanStimIm = stimImSum/(ceil(timePt5sec)-ceil(timePt3sec));
    stimMap = (meanStimIm-meanBaseIm)./meanBaseIm;
    filtBase = imgaussfilt(meanBaseIm, 2);
    filtStim = imgaussfilt(meanStimIm, 2);
    filtStimMap = (filtStim-filtBase)./filtBase;
    
    heatmap_fig = figure; colormap('jet');
    subplot(1,2,1); imagesc(stimMap); colorbar
    title('Avg dF/F0 During Stim (Raw)')
    subplot(1,2,2); imagesc(filtStimMap); colorbar
    title('Avg dF/F0 During Stim (Filtered)')
    
    tempDir = pwd;
    heatmap_figName = ['heatmap_' tempDir(strfind(tempDir,'TSer'):end)];
    try    
        saveas(heatmap_fig,heatmap_figName,'fig');
    catch
        disp('Problem with heatmap figure saving.'); 
    end
    
    %%% Save data to a file containing the Tseries name
    currentDirectory = pwd;
    TseriesLoc = strfind(currentDirectory, 'TSeries-');
    if isempty(TseriesLoc)
        TseriesLoc = strfind(currentDirectory, 'Tseries-');
    end
    cd ..
    parentDir = pwd;
    cd(currentDirectory); % Change back to the tseries-containing directory.
    
    saveName = NaN;
    while isnan(saveName)
        for end_idx = 1:20 %length(parentDir)
            if strcmp(parentDir(end-(end_idx)), '/') || strcmp(...
                    parentDir(end-(end_idx)), '\')
                saveName = parentDir(end-end_idx+1:end);
            end
        end
    end
    
    % Get trial type:
    if strcmp(parentDir(end-2:end), 'INE')
        trialType = strcat('_', 'PRE-EXPT_BASELINE');
    elseif strcmp(parentDir(end-2:end), 'T5V')
        trialType = strcat('_', 'PT5V');
    elseif strcmp(parentDir(end-3:end), '15')
        trialType = strcat('_', '15V');
    else
        trialType = strcat('_', parentDir(end-2:end));
    end
    
    timePoints_FileName = strcat(currentDirectory(TseriesLoc:end), trialType, '_timePoints');
    
    delFoverF_FileName = strcat(currentDirectory(TseriesLoc:end), trialType, '_delFoverF');
    
    baselineF0_FileName = strcat(currentDirectory(TseriesLoc:end), trialType, '_baselineF0');
    
    
    %% Save copy of each variable in the specified folder:
%     savePath = currentDirectory(...
%         strfind(currentDirectory, 'Raw_Un')+16:strfind(...
%         currentDirectory, 'eries')-4);
%     saveDirectory = strcat(prefixDirectory, '/Data_Analysis/', ...
%         savePath);  

    if useROI == 1 % if using ROI setting
        disp('ROI name is: '); disp(roiName);
        
        roiSaveFolderPre = strcat(exptFolder,'/',roiName,'/'); % to add a subfolder so that the different ROI analysis sets end up in their own subfolders
        roiSaveFolder = strcat(roiSaveFolderPre,stimFolder,'_',roiName);
        saveDirectory = strrep(roiSaveFolder,'Raw_Unprocessed','Data_Analysis');
        
        disp('saveDirectory is: '); disp(saveDirectory);
    else
        saveDirectory = strcat(stimFolder,exptFolder);
    end
    
    if ~isdir(saveDirectory)
        mkdir(saveDirectory);
    end
    cd(saveDirectory);
    
    save(timePoints_FileName, 'timePoints');
    save(delFoverF_FileName, 'delFoverF');
    save(baselineF0_FileName, 'F_0'); % Save the baseline (F_0)!
    
    if usePID == 1
        saveas(PID_fig,PID_figName,'fig');
    end
    
    cd(currentDirectory);
    
    % % % Previously: Save a copy of each variable in the containing folder:
    % cd(parentDir)
    % save(timePoints_FileName, 'timePoints');
    % save(delFoverF_FileName, 'delFoverF');
    % cd(currentDirectory)
    
    %%% Plot GCaMP average versus time.
    % repetitionPeriod = 0.818748; % Frame period : seconds/frame.
    % repetitionPeriod fetched automatically from config file, stored as framePerVal.
    % timePoints = [ 1:maxTime ] * repetitionPeriod;
    
    figure;
    % plot(timePoints, CaAvg);
    plot(timePoints, delFoverF);
    hold on
    title('deltaF/F vs Time')
    % % Add a line at timepoint = 3 sec
    yvals = get(gca, 'ylim');
    xvalsOn = [3,3];
    xvalsOff = [5,5]; % use 4 for Pipette stim, 5 for Odor Deliv.
    lineOn = plot(xvalsOn, yvals);
    lineOff = plot(xvalsOff, yvals);
    % Make the new vertical line at 3s a dash-dotted black line.
    set(lineOn,'Color','k','LineWidth', 1, 'LineStyle','-.')
    set(lineOff,'Color','k','LineWidth', 1, 'LineStyle','-.')
    hold off
    
    % % Annotate the point (3, delFoverF(3))
    % text(3, delFoverF(3),'\leftarrow (delFoverF(3sec))',...
    %      'HorizontalAlignment','left')
    % title('GCaMP Average vs Time')
    % % Add a line at timepoint = 3 sec
    % % ylimits = ylim;
    % % maxY = max(delFoverF);
    % % line([3 0.01],[3 maxY],[1 1],'Marker','.','LineStyle','-')
    
    displayLine = ['Analysis for  ' ...
        currentDirectory(strfind(currentDirectory, '201'):end) ...
        ' completed.'];
    disp(displayLine);
    
    %%% Demonstration code lines from Daisuke:
    %
    % ans =
    %
    % /Users/katherineshakman/Documents/MATLAB:
    %
    % tif = dir('*.tif');
    % im1 = imread(tif(1).name));
    %  im1 = imread(tif(1).name));
    %                          |
    % Error: Unbalanced or unexpected parenthesis or bracket.
    %
    % im1 = imread(tif(1).name);
    % figure;imshow(uint8(im1))
    % figure;imshow(uint8(im1./8))
    % x = impoly
    %
    % x =
    %
    %   impoly handle with no properties.
    %   Methods, Events, Superclasses
    %
    % pos = getPosition(x);
    % pos
    %
    % pos =
    %
    %    73.0000  355.0000
    %   146.0000  333.0000
    %   153.0000  431.0000
    %    91.0000  445.0000
    %
    % mask = createMask(x);
    % figure;imshow(mask)
    % roiv1 = im1.*mask;
    % Error using  .*
    % Integers can only be combined with integers of the same class, or scalar doubles.
    %
    % roiv1 = double(im1).*mask;
    % figure;imshow(uint8(roiv1./8))
    % figure;imshow(uint8(roiv1./12))
    % figure;imshow(uint8(roiv1./16))
    cd ../
    display('new_FluorAnalysis_batch completed');
end

