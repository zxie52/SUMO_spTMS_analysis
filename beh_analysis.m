
fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120', ...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};
       
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

group_name = {'left AMI', 'left UMI', 'right AMI', 'right UMI'};

% create emoty accuracy table for later data appending
accuracy = [];
accuracyProbe1 = [];
accuracyProbe2Stay = [];
accuracyProbe2Switch = [];

for l = 1:length(subject)
    % we will just use the cue1 eeg file to extrac the beh data on accuracy and response time
    for t = 2
        % load the EEG file
        EEG = pop_loadset('filename', strcat(subject{l}, '_before_iem_', type{t}, '.set'), 'filepath', fpath8);
        
        % change the EEG.event to a cell array and then categorize into different bins
        p = struct2table(EEG.event); 
        
        %% Calculate the mean for different conditions
        % times 100 to change mean to percentage
        % have the response on probe1
        probe1Accuracy = mean(p.response, 1) * 100;

        % have the response on probe2
        probe2AccuracyStay = mean(p.response2(p.targetlocation == p.targetlocation2),1) * 100;
        probe2AccuracySwitch = mean(p.response2(p.targetlocation ~= p.targetlocation2),1) * 100;

        % have the response for specific conditions in probe1(left item and right item)
        probe1AccuracyLeft = mean(p.response(p.targetlocation == 1), 1) * 100;
        probe1AccuracyRight = mean(p.response(p.targetlocation == 2), 1) * 100;
        
        % have the response for specific conditions in probe2 stay(left item and right item)
        % targetlocation and targetlocation2 are same and targetlocation2 is 1 or 2
        probe2AccuracyStayLeft = mean(p.response2(p.targetlocation == p.targetlocation2 & p.targetlocation2 == 1),1) * 100;
        probe2AccuracyStayRight = mean(p.response2(p.targetlocation == p.targetlocation2 & p.targetlocation2 == 2),1) * 100;

        % have the response for specific conditions in probe2 switch(left item and right item)
        probe2AccuracySwitchLeft = mean(p.response2(p.targetlocation ~= p.targetlocation2 & p.targetlocation2 == 1),1) * 100;
        probe2AccuracySwitchRight = mean(p.response2(p.targetlocation ~= p.targetlocation2 & p.targetlocation2 == 2),1) * 100;        
        
        %% Create the table for the accuracy values
        % create an accuracy table for probe1Accuracy, probe2AccuracyStay, and probe2AccuracySwitch
        accuracy = [accuracy; probe1Accuracy,probe2AccuracyStay, probe2AccuracySwitch];
        % create an accuracy table for probe1Accuracy, probe1AccuracyLeft and probe1AccuracyRight
        accuracyProbe1 = [accuracyProbe1; probe1Accuracy, probe1AccuracyLeft, probe1AccuracyRight];
        % create an accuracy table for probe2AccuracyStay, probe2AccuracyStayLeft and probe2AccuracyStayRight
        accuracyProbe2Stay = [accuracyProbe2Stay; probe2AccuracyStay, probe2AccuracyStayLeft, probe2AccuracyStayRight];
        % create an accuracy table for probe2AccuracyStay, probe2AccuracyStayLeft and probe2AccuracyStayRight
        accuracyProbe2Switch = [accuracyProbe2Switch; probe2AccuracySwitch, probe2AccuracySwitchLeft, probe2AccuracySwitchRight];
    end
end