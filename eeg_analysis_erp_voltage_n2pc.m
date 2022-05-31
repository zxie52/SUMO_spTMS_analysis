%% Load the Basics :: For Plotting the ERPs
%%%%%%%%%%%%%%%%%%%%% For ERP plotting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
addpath(genpath('E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\analysis_code'));
fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';
fpath14 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\ERP_results_voltage';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120', ...
           'SUMO_3001', 'SUMO_3017'};
       
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

group_name = {'left trials', 'right trials'};

for t = 1 : 7 % for seven intervals
    fprintf(strcat("We are running the ", type{t}, " interval for ERP analysis\n"));
    pause_script = input('Press enter to continue');
    for s = 1 : length(subject)
        %% Step 1: only pick up the correct trials for plotting the ERP
        cd(fpath8);
        EEG = pop_loadset('filename', strcat(subject{s}, '_before_iem_', type{t}, '.set'), 'filepath', fpath8);

        p = struct2table(EEG.event); % change the EEG.event to a cell array and then categorize into different bins
        % filter only correct trials in the probe1
        q = p((p.response == 1), :);

        %% Step 2: filter trials into different conditions
        if ismember(t, 1:4), targetloc = q.targetlocation;
        elseif ismember(t, 5:7), targetloc = q.targetlocation2;
        end
        
        left_trials = q.epoch(targetloc == 1);
        right_trials = q.epoch(targetloc == 2);
        groups = {left_trials, right_trials};
        
        %% Extra step before plotting: remove the baseline from -200ms to 0ms
        EEG = pop_rmbase( EEG, [-200 0]);
        
        %% Step 3: save out the eeg data in different conditions
        for i = 1:length(groups)
            %for contraposterior electrodes
            switch group_name{i}
                case {'left trials'}
                    ROI = {'CP2', 'CP4', 'CP6', 'TP8', 'TP10', ...
                           'P2', 'P4', 'P6', 'P8', ...
                           'PO4', 'PO8', 'O2'};
                    contra_impchan = find(ismember({EEG.chanlocs.labels}, ROI)); %channels in R hem
                    ROI = {'CP1', 'CP3', 'CP5', 'TP7', 'TP9', ...
                           'P1', 'P3', 'P5', 'P7', ...
                           'PO3', 'PO7', 'O1'};
                    ipsi_impchan = find(ismember({EEG.chanlocs.labels}, ROI)); %channels in L hem
                case {'right trials'}
                    ROI = {'CP1', 'CP3', 'CP5', 'TP7', 'TP9', ...
                           'P1', 'P3', 'P5', 'P7', ...
                           'PO3', 'PO7', 'O1'};
                    contra_impchan = find(ismember({EEG.chanlocs.labels}, ROI)); %channels in L hem
                    ROI = {'CP2', 'CP4', 'CP6', 'TP8', 'TP10', ...
                           'P2', 'P4', 'P6', 'P8', ...
                           'PO4', 'PO8', 'O2'};
                    ipsi_impchan = find(ismember({EEG.chanlocs.labels}, ROI)); %channels in R hem

            end
            % average the eeg data across trials then across impchans(contraposterior channels)
            eeg_data_contra = squeeze(mean(mean(EEG.data(contra_impchan, :, groups{i}),3),1));
            eeg_data_ipsi = squeeze(mean(mean(EEG.data(ipsi_impchan, :, groups{i}),3),1));
            % save the average in 2 condtions 
            subject_condition_contra(i,:) = eeg_data_contra;
            subject_condition_ipsi(i,:) = eeg_data_ipsi;
        end

        %% Step 4: Save data into one bigger matrix: subjects * timeseries * condtions
        eeg_group_contra(s, :, :) = subject_condition_contra;
        eeg_group_ipsi(s, :, :) = subject_condition_ipsi;
        dstimes = EEG.times;
        clear subject_condition;
    end

    %% Step 5: Plot the ERP figure for this condition
    % filter the data in gaussian filter for 16ms for plotting
    for g = 1 : size(eeg_group_contra, 2)
        for i = 1:size(eeg_group_contra, 1)
            for j = 1:size(eeg_group_contra,2)
                smoothed_contra(i,j,:) = smoothdata(eeg_group_contra(i,j,:), 'gaussian', 16);
            end
        end

        for i = 1:size(eeg_group_ipsi, 1)
            for j = 1:size(eeg_group_ipsi,2)
                smoothed_ipsi(i,j,:) = smoothdata(eeg_group_ipsi(i,j,:), 'gaussian', 16);
            end
        end

        % set up the figures for the ERP plotting
        [titleFontSize, axisFontSize, textFontSize] = erp_plotting_setup();

        % calculate the data that will be plotted
        data_ave_contra = squeeze(mean(smoothed_contra, 1));
        data_err_contra = squeeze(std(smoothed_contra,1)) / sqrt(size(smoothed_contra,1));

        % calculate the data that will be plotted
        data_ave_ipsi = squeeze(mean(smoothed_ipsi, 1));
        data_err_ipsi = squeeze(std(smoothed_ipsi,1)) / sqrt(size(smoothed_ipsi,1));

%         % remove the baseline before plotting(from -200ms to 0ms)
%         data_ave_contra = data_ave_contra - mean(data_ave_contra(:, find(dstimes == -200):find(dstimes == 0)), 2);
%         data_err_contra = data_err_contra - mean(data_err_contra(:, find(dstimes == -200):find(dstimes == 0)), 2);
%         data_ave_ipsi = data_ave_ipsi - mean(data_ave_ipsi(:, find(dstimes == -200):find(dstimes == 0)), 2);
%         data_err_ipsi = data_err_ipsi - mean(data_err_ipsi(:, find(dstimes == -200):find(dstimes == 0)), 2);
      
        % plot the trigger line in the figure
        ylim([-4, 4]);
        y = ylim;
        p = plot([0 0], [y(1) y(2)]);
        p.Color = 'black';
        p.LineWidth = 2;
        % Add mark onto the TMS trigger line
        txt = text(0, y(1)-.5, type{t}, 'Fontsize', textFontSize);
        txt.FontWeight = 'bold';
        txt.HorizontalAlignment = 'center';

        % plot the ERPs for four conditions in one figure
        % for left trials contra channels
        p1 = plot(dstimes, data_ave_contra(g,:), 'lineWidth', 3, 'Color', 'r');
        shadedError(dstimes, data_ave_contra(g,:), data_err_contra(g,:), 'r');
        % for left trials ipsi channels
        p2 = plot(dstimes, data_ave_ipsi(g,:), 'lineWidth', 3, 'Color', 'b');
        shadedError(dstimes, data_ave_ipsi(g,:), data_err_ipsi(g,:), 'b');

        % show the legend in the figure
        condition_name = {'contra', 'ipsi'};
        legend([p1, p2], condition_name, 'Fontsize', textFontSize);

        % write down the title, xlabel, and ylabel
        xlabel("Time across the trial Interval(ms)", 'Fontsize', axisFontSize);
        ylabel("Amplitude (µV)", 'Fontsize', axisFontSize);
        title(strcat("ERP on in ", type{t}, " interval(", group_name{g}, ")"), 'Fontsize', titleFontSize);
        cd(fpath14);
        savefig(strcat("ERP on in ", type{t}, " interval(", group_name{g}, ")"));
        saveas(gcf, strcat("ERP on in ", type{t}, " interval(", group_name{g}, ").png"));
    end
end