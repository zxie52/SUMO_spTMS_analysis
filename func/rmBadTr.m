function EEGupdated = rmBadTr(EEG, type, path, subject)    
% Basic idea is that we will preprocess the stim, cue1 and cue2, if the trials in stim, cue1 and cue2 that
% subject blink during the onset of the stimulus, then we will have to remove that trial for later all analysis

% Input: EEG dataset, type of the condtion, directory pathway to save the EEG.BadTrManual, and the subject ID
% Output: EEG data in which the bad trials have been excluded

nEpochs = length(EEG.epoch);
EEGupdated = EEG;
% Bad trials which were dectected by automatic algorithm
EEGupdated.BadTrAuto = unique(find(EEG.reject.rejjp==1));
% Bad trials which were dectected manually
EEGupdated.BadTrManual = unique(find(EEG.reject.rejmanual==1));

% create the empty matrix for storing the trial number for bad trials
badTrialManual = [];
cd(path);
switch type
    % if the subject blinks during the onset of the stimuls, then we need to exclude the whole trial
    % if the subject blinks during the onset of the cue1, then we need to exclude this trial on
    % cue1/tms1/probe1
    % if the subject blinks during the onset of the cue2, then we need to exclude this trial on
    % cue2/tms2/probe2
    % if the subject blinks during the onset of the probe1/probe2, then we need to exclude this
    % trial on probe1 or probe2
    
    case {'stim'}
        % Bad trials which were detected by manual inspection (e.g blinking during the onset of the cue)
        % if there are manually rejected trials, reject them and save the # of the trials
        if ~isempty(EEGupdated.reject.rejmanual)
            EEGupdated = pop_rejepoch( EEGupdated, [EEGupdated.BadTrManual EEGupdated.BadTrAuto] ,0);
            % save the bad trials that manually selected because of blinking during the onset of cue
            badTrialManual = EEGupdated.BadTrManual;
            save(strcat(subject, "_bad_trial_stim.mat"), 'badTrialManual');
        else
            EEGupdated = pop_rejepoch( EEGupdated, EEGupdated.BadTrAuto, 0);
        end
    case {'cue1'}
        % Bad trials which were detected by manual inspection (e.g blinking during the onset of the cue)
        % if there are rejected trials in the stimulus epoch, then load it
        if isfile(strcat(subject, "_bad_trial_stim.mat"))
            f = load(strcat(subject, "_bad_trial_stim.mat"));
        end
        if exist('f', 'var')
            badTrialManual = [f.badTrialManual, EEGupdated.BadTrManual];
        else 
            badTrialManual = [badTrialManual, EEGupdated.BadTrManual];
        end
        EEGupdated = pop_rejepoch( EEGupdated, [badTrialManual EEGupdated.BadTrAuto] ,0);
        
        % if there are trials that were rejected because of blinking
        if ~isempty(badTrialManual)
            save(strcat(subject, "_bad_trial_epoch1.mat"), 'badTrialManual');
        end
    case {'tms1', 'probe1'}
        % if there are trials rejected in stim or cue1
        if isfile(strcat(subject, "_bad_trial_epoch1.mat"))
            f = load(strcat(subject, "_bad_trial_epoch1.mat"));
        end
        % reject both manual rejected trials and auto rejected trials
        EEGupdated = pop_rejepoch( EEGupdated, [f.badTrialManual, EEGupdated.BadTrAuto], 0);
    case {'cue2'}
        % Bad trials which were detected by manual inspection (e.g blinking during the onset of the cue)
        % if there are rejected trials in the stimulus epoch, then load it
        if isfile(strcat(subject, "_bad_trial_stim.mat"))
            f = load(strcat(subject, "_bad_trial_stim.mat"));
            badTrialManual = [f.badTrialManual, EEGupdated.BadTrManual];
        else 
            badTrialManual = [badTrialManual, EEGupdated.BadTrManual];
        end
        % reject both manual rejected trials and auto rejected trials
        EEGupdated = pop_rejepoch( EEGupdated, [badTrialManual EEGupdated.BadTrAuto] ,0);
        % if there are trials that were rejected because of blinking
        if ~isempty(badTrialManual)
            save(strcat(subject, "_bad_trial_epoch2.mat"), 'badTrialManual');
        end
    case {'tms2', 'probe2'}
        % if there are trials rejected in stim or cue2
        if isfile(strcat(subject, "_bad_trial_epoch2.mat"))
            f = load(strcat(subject, "_bad_trial_epoch2.mat"));
        end
        % reject both manual rejected trials and auto rejected trials
        EEGupdated = pop_rejepoch( EEGupdated, [f.badTrialManual, EEGupdated.BadTrAuto], 0);
    otherwise
        error("The condtion is not included in this function");
end

fprintf("There are %d trials removed in all %d trials in subject: %s\n", length([badTrialManual, EEGupdated.BadTrAuto]), nEpochs, subject);
pause_script = input('To keep going, press enter!\n');





