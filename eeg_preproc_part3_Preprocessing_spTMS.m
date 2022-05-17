close all;
clear;

%% Start the Journey

% Generally, there will be a few sessions for the preprocessing:
% Step 1: Merge the beh data on to the raw eeg file(done on the CRC due to the size of raw eeg file)
% Step 2: Epoching the data file(done on the CRC due to the size of the raw eeg file)
% Step 3: Preprocessing the eeg files for the IEM anaylsis and ERP/ERSP analysis
% 
% After the preprocessing, we will go to the analysis
% This part is only for the TMS part of the EEG
% for other parts(i.e stim, cue, and probe, please use xxx_noTMS.m in this folder)
%
% This script requires EEGLAB v2020.0 in MATLAB

% The flow for this part:
%                               Loading EEG sets
%                                       |
%                    Remove pre-stimulus baseline to remove the DC offset
%                                       |
%                              Remove poor channels
%                                       |
%          Remove TMS pulse artifact and peaks fo TMS-evoked muscle activity(-15 to 20ms) 
%                                       |
%                       Interpolate missing data around TMS pulse
%                                       |
%                     Downsample data (10,000 Hz to 1,000 Hz)
%                                       |
%            Remove TMS pulse artifact and peaks fo TMS-evoked muscle activity(-15 to 20ms)
%                                       |
%                    Detrend the EEG.data for removing the remainng DC offset
%                                       |
%        High pass through 1 Hz before 1st round of ICA to remove the remaining DC offset
%                                       |
%             1st round of ICA to remove the remaining TMS-evoked muscle actitivy
%                                       |
%                      Bandpass Filter(1-100Hz, excluding 55-65Hz)
%                                       |
%               2nd round of ICA to remove other artifacts(e.g blinking)
%                                       |
%                Interpolate the missing data and removed poor channels
%                                       |
%                               Re-reference to median
%                                       |
%                                    Save out
%                                       
% Be careful: This script is only for preprocessing the TMS 1 or TMS 2 intervals
% Also, because the ICA needs interaction between frontend and user, it is highly recommended to run
% this preprocessing scipt locally

%% Loading basics

% before the start of the preprocessing, download the eeg files from CRC to local machine
fpath3 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_preproc';
fpath4 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_1stICA';
fpath7 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_2ndICA';
fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120',...
           'SUMO_3001', 'SUMO_3017'};
label = {'tms1', 'tms2'};

% trigger label from EEG vmrk files
trigger = {'T  1', 'T  2'};
       
% choose the condition, tms1 or tms2
type  = input("What type of condition do you want to preprocess? tms1 or tms2 \n", 's');
while ~strcmp(type, {'tms1', 'tms2'})
    warning("For non-tms intervals, please use the script xxx_noTMS.m");
    fprintf("Please type in the valid type: tms1 or tms2? \n");
    type  = input("What type of condition do you want to preprocess? tms1 or tms2 \n", 's');
end
%% Saving baseline for later tesa
for l = 1% : length(subject)
    %% Step 1: Load eeg files which were epoched around TMS 1
    cd(fpath3);
    EEG = pop_loadset('filename', strcat(subject{l}, '_after_epoch_', type, '.set'), 'filepath', fpath3);
    eeglab redraw;
    fprintf(strcat("Preprocessing for Subject: ", subject(l), ", has been started! \n"));

    %% Extra step: Remove extra labels in the EEG.event
    % remove the extra labels S 50 and S 53
    if strcmp(type, 'tms1')
        EEG.event(strcmp({EEG.event.type}, 'S 50')) = [];
    elseif strcmp(type, 'tms2')
        EEG.event(strcmp({EEG.event.type}, 'S 53')) = [];
    end
    eeglab redraw;
    
    %% Step 1: Remove the DC offset in -250ms to -20ms in the pre-TMS period
    % Since we do not run the TMS-picker for checking the TMS onset
    % this time
    %%% This part is only for TMS EEG PART
    EEG = pop_rmbase( EEG, [-250 -20]);
    % ^TESA Notes: demeaning data is an alternative to remove DC-offsets (since baseline correction can reduce ICA reliability)
    %MLW20180214: need to epoch the data before demeaning
    %% Extra Step: Possibly, there will be bad channels during recording 
    % save channel locations for later interpolation 
    EEG = pop_select( EEG,'nochannel',{'HEOG' 'VEOG' 'Photodiode'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    EEG.allchan = EEG.chanlocs;
    if strcmp(subject{l}, 'SUMO_0102')
        % has huge TMS artifcat on FT9
        EEG = pop_select( EEG,'nochannel',{'TP9'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    end
    if strcmp(subject{l}, 'SUMO_0105')
        % has huge TMS artifcat on FT9
        EEG = pop_select( EEG,'nochannel',{'P3'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    end
    if strcmp(subject{l}, 'SUMO_0108')
        % has huge TMS artifcat on FT9
        EEG = pop_select( EEG,'nochannel',{'TP9'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    end
    if strcmp(subject{l}, 'SUMO_0120')
        % has huge TMS artifcat on FT9
        EEG = pop_select( EEG,'nochannel',{'FT9'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    end
    if strcmp(subject{l}, 'SUMO_3001')
        % has huge TMS artifcat on F1 and F3
        EEG = pop_select( EEG,'nochannel',{'F3' 'F1'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    end

    %% Step 2: Remove TMS pulse artifact and peaks fo TMS-evoked muscle activity(-15 to 20ms) -- TESA
    % Due to the fact that the TMS trigger is off for many trials, we remove the EEG data
    % between 15ms to 20ms to remove most of the TMS activity.
    EEG = pop_tesa_removedata( EEG, [-15 20] );

    %% Step 3: Interpolate missing data around TMS pulse - from TESA
    EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );
    %^TESA Notes: cubic interoplation replaces the removed data - needed to
    %help minimize ringing/step artifacts caused by the low-pass filter applied prior to downsampling

    %% Step 4: Downsample data (10,000 Hz to 1,000 Hz) - from TESA -- need to be 1000Hz so sampling rates match for TEP extraction
    EEG = pop_resample(EEG, 1000);
    eeglab redraw;

    %% Step 5: Replace interpolated data around TMS pulse with constant amplitude data (-3 to 10 ms) - from  TESA
    EEG = pop_tesa_removedata( EEG, [-15 20] );

    % Detrend the EEG.data for removing the DC offset
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end
    eeglab redraw

    % High pass through 1 Hz before 1st round of ICA
    EEG = pop_eegfiltnew(EEG, 'locutoff',1);
    eeglab redraw;

    %% Extra Step:  Remove bad trials - from TESA (run whole chunk at a time)
    EEG = rmBadTr(EEG, type, fpath8, subject{l});
    
    % save out before the 1st round of ICA
    cd(fpath4);
    EEG = pop_saveset(EEG, strcat(subject{l}, '_before_1stICA_', type, '.set'));

    %% Step 6: 1st ICA
    % The purpose for the 1st ICA is to remove the TMS-evoked muscle activity

    % stabilization is off for spTMS
    % stabilization is on for ppTMS(SICI, ICF, LICI)
    stabilization = 'off';
    EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' ); %set to on for ppTMS, off for spTMS
    EEG = pop_tesa_compselect( EEG,'comps',40,'figSize','medium','plotTimeX',[-100 400],'plotFreqX',[1 100],'tmsMuscle','on','tmsMuscleThresh',8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off','blink','off','blinkThresh',2.5,'blinkElecs',{'Fp1', 'Fp2'},'blinkFeedback','off','move','off','moveThresh',2,'moveElecs',{'F7', 'F8'},'moveFeedback','off','muscle','off','muscleThresh',0.6,'muscleFreqWin',[30 100],'muscleFeedback','off','elecNoise','off','elecNoiseThresh',4,'elecNoiseFeedback','off' );
    % In the ICA results, the TMS-evoked muscle may not be the first component
    % The DC offsets are not removed according to the ICA result
    % Perhaps we also need to add detrend process before the ICA
    fprintf(strcat(subject{l}, '\n'));

    %% Step 7: Interpolate missing data around TMS pulse - from  TESA
    EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );

    %% Step 8: Bandpass (1-100 Hz) and bandstop (55-65 Hz) filter data - from  TESA
    EEG = pop_tesa_filtbutter( EEG, 1, 100, 4, 'bandpass' ); %bandpass 1-50Hz like in Julkunen paper
    EEG = pop_tesa_filtbutter( EEG, 55, 65, 4, 'bandstop' ); %bandpass 1-50Hz like in Julkunen paper
    % In the SUMO, we only picked the bandpass (1-50Hz), some people used 1-80Hz

    %% Step 9: Replace interpolated data around TMS pulse with constant amplitude data
    EEG = pop_tesa_removedata( EEG, [-15 20] );

    % Remove first 100ms and last 100ms to avoid the border effect after filtering
    EEG = pop_epoch( EEG, trigger(strcmp(type, label)), [-0.2 .8]);
    eeglab redraw;
    
    % save out the eeg before the 2nd round of ICA
    cd(fpath7);
    eeglab redraw;
    EEG = pop_saveset(EEG, strcat(subject{l}, '_before_2ndICA_', type, '.set')); 

    %% Step 10: 2nd ICA
    % The purpose of 2nd ICA is same with normal ERP preprocessing, to remove artifcats, especially eye
    % blinks(even though we do not really care bout the frontal channels)

    % for SUMO102-120, eye blinks were picked by Fp1 and Fp2; eye movements were picked by F7&F8
    EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
    EEG = pop_tesa_compselect( EEG,'compCheck','on','comps',40,'figSize','medium','plotTimeX',[-100 400],'plotFreqX',[1 100],...
                               'tmsMuscle','on','tmsMuscleThresh',8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off',...
                               'blink','on','blinkThresh',2.5,'blinkElecs',{'Fp1' 'Fp2'},'blinkFeedback','off',...
                               'move','on','moveThresh',2,'moveElecs',{'Fp1' 'Fp2'},'moveFeedback','off',...
                               'muscle','on','muscleThresh',0.6,'muscleFreqIn',[7 75], 'muscleFreqEx',[58 62],'muscleFeedback','off',...
                               'elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );
    %% Step 11: Interpolate missing data around TMS pulse - from  TESA
    EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );

    % If there are channels removed before the preprocessing
    EEG = pop_interp(EEG, EEG.allchan, 'spherical');

    %% Step 12: Re-reference to median
    %median re-reference - better than average b/c average sensitive to outliers
    EEG.data=bsxfun(@minus,EEG.data,median(EEG.data,1)); 
    eeglab redraw

    %% Step 13: Save out the eeg dataset for later IEM analysis
    cd(fpath8);
    EEG = pop_saveset(EEG, strcat(subject{l}, '_before_iem_', type, '.set')); 
    fprintf(strcat("Preprocessing for Subject: ", subject(l), ", has been finished! \n"));
end 

fprintf(strcat("All preprocessings have been finished! \n"));
