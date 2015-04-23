% Make folders for imaging expts
BASIC_NAME = 'TRIAL';
NUM_TRIAL_TYPES = 6;
for idx = 60:70 %50
    
    %if idx > 5 && idx < 36
    %     if idx>47 && idx <60
    %         if mod(idx,NUM_TRIAL_TYPES) == 0
    %             trialType = '_PT5V';
    %         elseif mod(idx,NUM_TRIAL_TYPES) == 1
    %             trialType = '_10V';
    %         elseif mod(idx,NUM_TRIAL_TYPES) == 2
    %             trialType = '_1V';
    %         elseif mod(idx,NUM_TRIAL_TYPES) == 3
    %             trialType = '_8V';
    %         elseif mod(idx,NUM_TRIAL_TYPES) == 4
    %             trialType = '_2V';
    %         elseif mod(idx,NUM_TRIAL_TYPES) == 5
    %             trialType = '_4V';
    %         end
    %     else
    %         trialType = '_BASELINE';
    %     end
    if idx > 59
        trialType = '_BASELINE';
    end
    newFolderName = [BASIC_NAME, num2str(idx), trialType];
    mkdir(newFolderName);
    disp(newFolderName)
    disp('---')
end