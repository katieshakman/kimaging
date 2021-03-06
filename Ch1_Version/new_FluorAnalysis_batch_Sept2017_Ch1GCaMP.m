%Katie Shakman 
%% Fluorescence Analysis
% Code to process GCaMP calcium imaging Tseries (from 2-photon).
% What it does/will do:
% Reads image, allows selection of an ROI, gets fluorescence information
% for the ROI over time, saves a chart of the deltaF/F values, and makes a
% plot of deltaF/F for the ROI over time.

% Imaging channels for 9th Floor:
% Ch1 : green gaasp (GCaMP or GFP)
% Ch2 : red PMT 1 HV (RFP)

clear all; close all; 

%% User Parameters/Settings

getAllChAvg = 1; % Set to 1 to get an average image in all channels. 
getCh1Avg = 0; % Set this value to 1 to get the Ch1 Avg Img.
useROI = 1; % set useROI to 1 to load and apply an ROI  
stim = 'Ch1'; 
stimOn = 6; 
stimOff = 8; 
TSvsCycle = 0; 

startDir = pwd; % default starting dir

% Indicate whether to use long format directory name for saving dir.
LongSavePath = 1; 

%% Load ROI
if useROI == 1
    [ROIfile, ROIfilePath] = uigetfile('Please select ROI file'); 
    cd(ROIfilePath); 
    load(ROIfile); % loads as variable 'ROI'
    cd(startDir); 
    roiName = ROIfile(1:end-4);
end
%% Compile list of Tseries directories in current folder.
myFolder = pwd; 
if TSvsCycle == 0
    TS_listing = dir('Tseries*');
    if isempty(TS_listing)
        TS_listing = dir('TSeries*');
        if isempty(TS_listing)
            display('No TSeries found in directory.')
        end
    end
elseif TSvsCycle == 1
    TS_listing = dir('Cycle*');
    if isempty(TS_listing)
        TS_listing = dir('Cycle*');
        if isempty(TS_listing)
            display('No Cycles found in directory.')
        end
    end
end
%% Read the image (and get an ROI, though this not implemented now).
for TS_idx = 1:length(TS_listing)
    close all; 
    % clear all variables except myFolder, TS_listing, and TS_idx:
    %clearvars -except myFolder TS_listing TS_idx LongSavePath getAllChAvg getCh1Avg useROI startDir; 
%     clear all; close all;
%     imDir = uigetdir(); % Allow user to specify the directory to read images from.
try 
    mySubfolder = TS_listing(TS_idx).name
    imDir = fullfile(myFolder, mySubfolder); 
    cd(imDir);
catch
    mySubfolder = TS_listing(TS_idx).name
    imDir = fullfile(myFolder, mySubfolder); 
    cd(imDir);
end
    tifIm = dir('*.tif');
    %im1 = imread(tifIm(1).name); % Below: changed this to the average image for the stack.
    cd ../; % Move to containing folder (up one level).
    currFolderName = pwd; % Get name of containing folder.
    cd(imDir);
    
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
%         clear runningTotIm;
    end
    
    % Specify ROI if desired
    if useROI == 1
        if ~exist('ROI') || isempty('ROI')
            ROI = uint16(ones(size(avgIm))); % Use entire average image size for ROI if NaN ROI was passed.
            roiName = 'full';
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
    else
        ROI = ones(size(avgIm));
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
        avgIm = runningTotIm./(length(tifIm)/3); % TODO: Find the actual # of Ch1 Ims.
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
    
    %% Get GCaMP average at each time point; convert to delF/F.
    % TODO: Add in ROI restriction.
    % Do this only for the part of tifIm that is Ch1 (GCaMP for 9th floor settings)
    
    %Ch2tifIm = %all of tifIm for which tifIm(1).name contains Ch2
    % TODO: Check that there are Ch1, Ch2, and Ch3 images
    % before dividing by 3 (if only Ch1 and Ch2, divide by 2).
    % % % maxTime = length(tifIm)/3; % Last timepoint/slice,
    % div by 3 since only second 1/3 of the Tseries files will be Ch2 (if Ch3 data exists).
    % TODO: Change to Ch1tifIm instead of tifIm, so we are sure we only use Ch1
    % images.
    
    tifIm_names = cell(length(tifIm),1);
    for idx = 1:length(tifIm)
        tifIm_names{idx} = tifIm(idx).name;
    end
    Ch1_indices = strfind(tifIm_names, 'Ch1'); % get the GC channel images
    tifImCh1_indices = cell2mat(Ch1_indices);
    % tifImCh1_indices = tifImCh1_indices * ones(length(tifImCh1_indices));
    maxTime = length(tifImCh1_indices);
    CaAvgDenom = size(ROI,1); % Total number of pixels per slice (usually). Check actual resolution in xml file.
    % TODO: Change the above (CaAvgDenom) to get the true image (or ROI) size.
    
    for timeIdx = 1:maxTime
        currTif = imread(tifIm(timeIdx).name) ;
        % Mask all images in stack and avg image with ROI (if using):
        if useROI
            currTif = ROI.*currTif;
        end
        CaAvg(timeIdx) = sum(sum( currTif )) / CaAvgDenom; % Avg Signal in frame normalized to the image size.
    end
    timePtStimOnsec = stimOn/framePerVal; % 3sec*(1/(sec/frame)) = 3sec/framePerVal
    if timePtStimOnsec < 1 || isnan(framePerVal)
        F_0 = mean(CaAvg(1:1)); % Average fluorescence in first frame.
    else
        F_0 = mean(CaAvg(1:floor(timePtStimOnsec))); % Average fluorescence up to 3 sec (takes floor of timePt3sec).
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
    if TSvsCycle == 0 
        TseriesLoc = strfind(currentDirectory, 'TSeries-');
        if isempty(TseriesLoc)
            TseriesLoc = strfind(currentDirectory, 'Tseries-');
        end
    elseif TSvsCycle == 1
        TseriesLoc = strfind(currentDirectory, 'Cycle');
        if isempty(TseriesLoc)
            TseriesLoc = strfind(currentDirectory, 'Cycle');
        end
    end
    cd ..
    parentDir = pwd;
    cd(currentDirectory); % Change back to the tseries-containing directory.
    
    saveName = NaN;
    tic 
    while isnan(saveName)
        for end_idx = 1:50 %length(parentDir)
            if strcmp(parentDir(end-(end_idx)), '/') || strcmp(...
                    parentDir(end-(end_idx)), '\')
                saveName = parentDir(end-end_idx+1:end);
            end
        end
        loopTime = toc;
        if loopTime > 10
            saveName = strcat(parentDir,'saveNameError'); 
            break
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
        stimFolder = stim; 
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
    
    
    try    
        save(timePoints_FileName, 'timePoints');
    catch
        save('timePoints','timePoints'); 
    end
    
    try
        save(delFoverF_FileName, 'delFoverF');
    catch
        save('delFoverF','delFoverF'); 
    end

    try
        save(baselineF0_FileName, 'F_0'); % Save the baseline (F_0)! 
    catch 
        save('F_0','F_0'); 
    end
    
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
    xvalsOn = [stimOn,stimOn]; 
    xvalsOff = [stimOff,stimOff]; % use 4 for Pipette stim, 5 for Odor Deliv.
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
% could also go back to the starting directory with: 
% cd(startDir); 