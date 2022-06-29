% generate the behavioral data to csv file
addpath(genpath('E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\analysis_code'));
fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';
fpath11 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_beh';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120',...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};
        
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

fprintf("You have following epochs can choose: %s, %s, %s, %s, %s, %s, %s\n", type{:});
condition = input("What epoch do you want to generate the beh data?\n", 's');
while ~ismember(condition, type)
    warning("The input epoch is not available, please re-type the epoch.");
    condition = input("Do you want the IEM on testorient(test orientation) or targetorient(target orientation)?\n", 's');
end

for s = 1 : length(subject)
    cd(fpath8);
    EEG = pop_loadset('filename', strcat(subject{s}, '_before_iem_', condition, '.set'), 'filepath', fpath8); 
    p = struct2table(EEG.event); % change the EEG.event to a cell array and then categorize into different bins
    
    cd(fpath11);
    writetable(p, strcat(subject{s}, "_beh_", condition, ".csv"));
end
