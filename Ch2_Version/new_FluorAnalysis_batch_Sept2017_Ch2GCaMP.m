%Katie Shakman 
%% Fluorescence Analysis
% Code to process GCaMP calcium imaging Tseries (from 2-photon).
% What it does/will do:
% Reads image, allows selection of an ROI, gets fluorescence information
% for the ROI over time, saves a chart of the deltaF/F values, and makes a
% plot of deltaF/F for the ROI over time.

% Imaging channels:
% Ch1 : red PMT 1 HV (RFP)
% Ch2 : green gaasp (GCaMP or GFP) on the 5th floor settings
% Ch3 : dodt (transmitted light)
clear all; close all; 

%% User Parameters/Settings

getAllChAvg = 1; % Set to 1 to get an average image in all channels. 
getCh2Avg = 0; % Set this value to 1 to get the Ch1 Avg Img.
useROI = 1; % set useROI to 1 to load and apply an ROI  
stim = 'Ch2'; 

startDir = pwd; % default starting dir

% Indicate whether to use long format directory name for saving dir.
LongSavePath = 1; 
numChannels = 2; 
imRes = 128*128; % could also read this from the XML file
endBaseline = 3; % assuming baseline calculation should end at time=3sec
%% Load ROI
if useROI == 1
%     ROIdir = uigetdir('Please select ROI directory'); 
    [ROIfile, ROIfilePath] = uigetfile('Please select ROI file'); 
    cd(ROIfilePath); 
    load(ROIfile); % loads as variable 'ROI'
    cd(startDir); 
    roiName = ROIfile(1:end-4);
end
%% Compile list of Tseries directories in current folder.
myFolder = pwd; 
TS_listing = dir('Tseries*');
if isempty(TS_listing)
    TS_listing = dir('TSeries*');
    if isempty(TS_listing)
        display('No TSeries found in directory.')
    end
end
%% Read the image (and get an ROI, though this not implemented now).
for TS_idx = 1:length(TS_listing)
    close all; 
%     imDir = uigetdir(); % Allow user to specify the directory to read images from.
    mySubfolder = TS_listing(TS_idx).name
    imDir = fullfile(myFolder, mySubfolder); 
    cd(imDir);
    tifIm = dir('*.tif');
    %im1 = imread(tifIm(1).name); % Below: changed this to the average image for the stack.
    cd ../; % Move to containing folder (up one level).
    currFolderName = pwd; % Get name of containing folder.
    cd(imDir);
    
    if getAllChAvg == 1
        % Get the average image for the stack in all channels:
        imSize = size(imread(tifIm(1).name));
        runningTotImAllCh = uint16(zeros(imSize)); % Initialize
        for idx = 1:length(tifIm)
            newTifIm = imread(tifIm(idx).name);
            lastIdx = idx; % Save the last idx value for the average.
            runningTotImAllCh = runningTotImAllCh + newTifIm;
        end
        % avgIm = zeros(512,512); % Initialize the average image variable.
        avgImAllCh = runningTotImAllCh./length(tifIm);
        %         clear runningTotIm;
    end
    
    % Specify ROI if desired (loads from file if the user selected one)
    if useROI == 1
        if ~exist('ROI') || isempty('ROI')
            ROI = uint16(ones(size(avgIm))); % Use entire average image size for ROI if NaN ROI was passed.
            roiName = 'full';
        elseif isnan(ROI);
            ROI = uint16(ones(size(avgIm))); % Use entire average image size for ROI if NaN ROI was passed.
            roiName = 'full';
        else
            figure;
            subplot(1,2,1); imshow(runningTotImAllCh);
            subplot(1,2,2); imshow(ROI.*(runningTotImAllCh));
            % TODO: Add a dialog to confirm or adjust ROI.
        end
        % Mask average image with ROI.
        avgImAllCh = ROI.*avgImAllCh;
    else
        ROI = ones(size(avgImAllCh));
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
    
    if getCh2Avg == 2
        % Get the average image for the stack in only Ch2:
        frameSize = imSize; % Auto-adjusts for size of frames in the Tseries.
        runningTotImCh2 = uint16(zeros(frameSize)); % Initialize
        pattern = 'Ch2'; % Pattern to search for in image name string.
        for idx = 1:length(tifIm)
            % Locate "Channel 2" in current image name:
            k = strfind(tifIm(idx).name, pattern);  % Returns empty var if pattern not found.
            if ~isempty(k) % Read im if k is NOT empty (Read Channel 2 images).
                newTifIm = imread(tifIm(idx).name);
            end
            runningTotImCh2 = runningTotImCh2 + newTifIm;
        end
        avgImCh2 = runningTotImCh2./(length(tifIm)/numChannels); % TODO: Find the actual # of Ch2 Ims.
        clear runningTotImCh2;
        clear k;
    end
    
    % Make the figure(s).
    % figure;
%     divImBy = 2; % Daisuke initially suggested 16 for this value.
    
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
    
    %% Get GCaMP average at each time point; convert to delF/F.
    % TODO: Add in ROI restriction.
    % Do this only for the part of tifIm that is Ch2 (GCaMP for 5th floor settings)
    
 
    % TODO: Check that there are Ch1, Ch2, and Ch3 images
    % before dividing by 3 (if only Ch1 and Ch2, divide by 2).
    % % % maxTime = length(tifIm)/3; % Last timepoint/slice,
    % div by 3 since only second 1/3 of the Tseries files will be Ch2 (if Ch3 data exists).
    
    tifIm_names = cell(length(tifIm),1);
    for idx = 1:length(tifIm)
        tifIm_names{idx} = tifIm(idx).name;
    end
    Ch2_indices = strfind(tifIm_names, 'Ch2');
    tifImCh2_counter = cell2mat(Ch2_indices);
    % tifImCh1_indices = tifImCh1_indices * ones(length(tifImCh1_indices));
    startAt = length(tifIm_names) - length(tifImCh2_counter) +1; % e.g. start at 1 if no Ch1, start at 101 if there are 100 frames of Ch1.
    maxTime = length(tifImCh2_counter);
    CaAvgDenom = imRes; % Total number of pixels per slice.
    % Can change the above (CaAvgDenom) to get the true image size
    % for an ROI if necessary: 
    if useROI
        CaAvgDenom = size(ROI,1)*size(ROI,2); 
    end
    
    CaAvg = []; % initialize. 
    for timeIdx = startAt:(startAt+maxTime-1) % e.g. go from frame 101 to 200.  Makes sure to skip the Ch1 frames if we are imaging with GAASP as Ch2. 
        currTif = imread(tifIm(timeIdx).name) ;
        % Mask all images in stack and avg image with ROI (if using):
        if useROI
            currTif = ROI.*currTif;
        end
        currCaAvg = sum(sum( currTif )) / CaAvgDenom; % Avg Signal in frame normalized to the image size.
        CaAvg = [CaAvg, currCaAvg]; 
    end
    timePtEndBaseline = endBaseline/framePerVal; % 3sec*(1/(sec/frame)) = 3sec/framePerVal, if baseline ends at 3sec
    if timePtEndBaseline < 1 || isnan(framePerVal)
        F_0 = mean(CaAvg(1:1)); % Average fluorescence in first frame.
    else
        F_0 = mean(CaAvg(1:timePtEndBaseline)); % Average fluorescence up to end of baseline (was 3 sec originally) 
    end
    
    % TODO: Vectorize below.
    for timeIdx = 1:maxTime 
        delF(timeIdx) = CaAvg(timeIdx) - F_0;
        delFoverF(timeIdx) = delF(timeIdx)/F_0;
    end
    timePoints = [ 1:maxTime ] .* framePerVal;
    % Above: time is in sec, framePerVal is in sec/frame --> divide by framePer
    % allTimePoints = [framePerVal:maxTime];
   
    %% Save data to a file containing the Tseries name
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
    
    if strcmp(parentDir(end-2:end), 'INE')
        trialType = strcat('_', 'PRE-EXPT_BASELINE');
    elseif strcmp(parentDir(end-2:end), 'T5V')
        trialType = strcat('_', 'PT5V');
    elseif strcmp(parentDir(end-3:end), '15')
        trialType = strcat('_', '15V');
    else
        trialType = strcat('_', parentDir(end-2:end));
    end
    % timePoints_FileName = strcat(currentDirectory(TseriesLoc:end), currFolderName, '_timePoints');
    
    % timePoints_FileName = [currentDirectory(TseriesLoc:end) trialType '_' 'timePoints'];
    timePoints_FileName = strcat(currentDirectory(TseriesLoc:end), trialType, '_timePoints');
    
    delFoverF_FileName = strcat(currentDirectory(TseriesLoc:end), trialType, '_delFoverF');
    
    baselineF0_FileName = strcat(currentDirectory(TseriesLoc:end), trialType, '_baselineF0');
    
    suffixDirectory = currentDirectory(...
        strfind(currentDirectory, 'Raw_Un')+16:...
        strfind(currentDirectory, 'lobe')+5);
    prefixDirectory = currentDirectory(...
        1:strfind(currentDirectory, 'Raw_Un')-2);
    % % New: Save copy of each variable in the user-specified folder:
    % saveDirectory = uigetfolder();
    if LongSavePath == 0
        saveDirectory = strcat(prefixDirectory, '/Data_Analysis/', ...
            suffixDirectory, saveName);
    else
        savePath = currentDirectory(...
            strfind(currentDirectory, 'Raw_Un')+16:strfind(...
            currentDirectory, 'eries')-4);
        saveDirectory = strcat(prefixDirectory, '/Data_Analysis/', ...
            savePath);
    end
    
    % Make new Data_Analysis subfolder to save everything into. 
    if useROI == 1 % if using ROI setting
        if ~exist('roiName')
            roiName = input('Please enter a name for the roi: ','s'); 
        end
        disp('ROI name is: '); disp(roiName);
        
        if ~exist('exptFolder')
            exptFolder = startDir;
        end
        roiSaveFolderPre = strcat(exptFolder,'/',roiName,'/'); % to add a subfolder so that the different ROI analysis sets end up in their own subfolders
        stimFolder =stim; 
        roiSaveFolder = strcat(roiSaveFolderPre,stimFolder,'_',roiName);
        saveDirectory = strrep(roiSaveFolder,'Raw_Unprocessed','Data_Analysis');
%         saveDirectory = strrep(roiSaveFolderPre,'Raw_Unprocessed','Data_Analysis');
          
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
    
    cd(currentDirectory);
    % % % Previously: Save a copy of each variable in the containing folder:
    % cd(parentDir)
    % save(timePoints_FileName, 'timePoints');
    % save(delFoverF_FileName, 'delFoverF');
    % cd(currentDirectory)
   
    %% Plot GCaMP average versus time.
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
end  
cd(saveDirectory); 
% cd('..'); % Up one level. 
% could also go back to the starting directory with: 
% cd(startDir); 