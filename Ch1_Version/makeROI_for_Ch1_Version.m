% makeROI

%% User parameters
chooseROI = 1; % default choice, asks for folder of .tiffs and loads ROI drawing interface. 
saveFilenameBase = 'ROI_'; 
%% Make the ROI
startDir = pwd; 
if chooseROI == 1
    [ROI, roiName] = drawROI();
else
    ROI = NaN; roiName = [];
end

%% Save the ROI
cd(startDir); % change back to the starting dir to save ROI there
saveFilename = [saveFilenameBase,roiName]; 
save(saveFilename,'ROI'); 