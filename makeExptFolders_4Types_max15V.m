% Make folders for imaging expts
BASIC_NAME = 'TRIAL';
NUM_TRIAL_TYPES = 4;
for idx = 1:35
    
    if idx > 5 && idx < 26
%       if idx < 26
            if mod(idx,NUM_TRIAL_TYPES) == 0
                trialType = '_15V';
            elseif mod(idx,NUM_TRIAL_TYPES) == 1
                trialType = '_2V';
            elseif mod(idx,NUM_TRIAL_TYPES) == 2
                trialType = '_8V';
            elseif mod(idx,NUM_TRIAL_TYPES) == 3
                trialType = '_4V';
            end
    else
        trialType = '_BASELINE';
%     if idx > 59
%         trialType = '_BASELINE';
%     end
    end
    newFolderName = [BASIC_NAME, num2str(idx), trialType];
    mkdir(newFolderName);
    disp(newFolderName)
    disp('---')
end