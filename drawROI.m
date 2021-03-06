function [ROI,roiName] = drawROI(max)

% Choose whether to display ch1, ch2, or all channels. 
% dispCh = 0; % to display all channels
dispCh = 1; % to display ch1 only
% dispCh = 2; % to display ch2 only 

origDir = pwd;  % Remember where we started. 
DispRange = [0, max]; % DisplayRange parameters to use (min, max)

% Get folder of images to use for ROI
dialog = 'Choose a folder of raw tif images on which to load for ROI selection:';
disp(dialog)
rawImsDir = uigetdir(origDir,dialog); 
cd(rawImsDir);  
% Select images from the appropriate channels to display. 
if dispCh == 1
    tifIm = dir('*Ch1*.tif');
elseif dispCh == 2
    tifIm = dir('*Ch2*.tif'); 
else
    tifIm = dir('*.tif');
end
%im1 = imread(tifIm(1).name); % Below: changed this to the average image for the stack.
%     cd ../; % Move to containing folder (up one level).
%     currFolderName = pwd; % Get name of containing folder.
%     cd(imDir);

getAllChAvg = 1;
    if getAllChAvg == 1
        % Get the average image for the stack in all specified channels:
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
    
    % Make ROI:
    fig = figure; 
    imShown = imshow(runningTotIm,'DisplayRange',DispRange);
    disp('Please draw an ROI on the figure and double-click on figure when done.') 
    ROI = roipoly;
    ROI = uint16(ROI);
    close;
    
    % Get name for ROI: 
    roiName = input('Input name for ROI: ','s');
    
    % Move back to original folder: 
    cd(origDir); 
