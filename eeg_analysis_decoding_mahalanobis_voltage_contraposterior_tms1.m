%% Load the Basics :: For TMS1 decoding
%%%%%%%%%%%%%%%%%%%%% For TMS1 decoding  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
addpath(genpath('E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\analysis_code'));
fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120', ...
            'SUMO_3017'};
% 'SUMO_3001',       
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

group_name = {'left AMI', 'left UMI', 'right AMI', 'right UMI'};


for l = 1:length(subject)
    for t = 3 % This time we only preprocess the Probe1 data
        %% Step 1: filter trials into four groups: left_AMI, left_UMI, right_AMI, right_UMI
        cd(fpath8);
        EEG = pop_loadset('filename', strcat(subject{l}, '_before_iem_', type{t}, '.set'), 'filepath', fpath8);

        p = struct2table(EEG.event); % change the EEG.event to a cell array and then categorize into different bins
        
        nbins = 7;
        binedges = round(linspace(1,181,nbins+1));
        bincent = round(mean([binedges(1:end-1);binedges(2:end)]));
        
        % filter only correct trials in the probe1
        q = p((p.response == 1), :);
        
        % filter trials into seven bins for later IEM
        [left_bin_AMI, left_bin_UMI, right_bin_AMI, right_bin_UMI] = filter_bins_for_iem_tms(q, binedges, bincent);

        %% Step 2: Preprocess data befor calculating Mahalanobis Distance
        groups = {left_bin_AMI, left_bin_UMI, right_bin_AMI, right_bin_UMI};
        for i = 1 : length(groups)
            % have the EEG data and stimlabels for the iem function
            % the output h: epoch# * stimlabels, only the right trials this time
            h = [];
            h = nonzeros(groups{i});
            [~, colIdcs] = find(groups{i} ~= 0);
            h(:,2) = bincent(colIdcs);
            h(h(:,1)>size(EEG.data, 3)/2, 2) = h(h(:,1)>size(EEG.data, 3)/2, 2) - 180;
            % sort the h to have the ordered stimlabels
            h = sortrows(h);

            %for contraposterior electrodes
            switch i
                case {1 2}
                    ROI = {'CP2', 'CP4', 'CP6', 'TP8', 'TP10', ...
                           'P2', 'P4', 'P6', 'P8', ...
                           'PO4', 'PO8', 'O2'};
                    impchan = find(ismember({EEG.chanlocs.labels}, ROI)); %channels in R hem
                otherwise
                    ROI = {'CP1', 'CP3', 'CP5', 'TP7', 'TP9', ...
                           'P1', 'P3', 'P5', 'P7', ...
                           'PO3', 'PO7', 'O1'};
                    impchan = find(ismember({EEG.chanlocs.labels}, ROI)); %channels in L hem
            end
            
            stimlabels = circ_ang2rad(h(:,2));
            super_charge = EEG.data(impchan,:,h(:,1));

            %% Calculate the Mahalanobis Distance
            data = single(permute(super_charge, [3 1 2]));
            theta = circ_ang2rad(h(:,2));
            n_folds = 2;

            % Ester
            [distance_cos,distances] = mahal_func_theta_kfold(data,theta,n_folds);
            distance_cos_group(l,:,i) = squeeze(mean(distance_cos,1));
        end
    end
end

%% Process the distance cosine in group
% calcualte the average and the SEM
distance_cos_group_ave = squeeze(mean(distance_cos_group, 1))';
distance_cos_group_err = squeeze(std(distance_cos_group,1) / sqrt(size(distance_cos_group,1)))';

% % smooth the data with Gaussian Filter time windowed in 16ms
% for i = 1:size(distance_cos_group_ave)
%     distance_cos_group_ave_smoothed(i,:) = smoothdata(distance_cos_group_ave(i,:), 'gaussian', 16);
%     distance_cos_group_err_smoothed(i,:) = smoothdata(distance_cos_group_err(i,:), 'gaussian', 16);
% end
% group_name = {'left AMI', 'left UMI', 'right AMI', 'right UMI'};
% plot the decoding results in the same figure
figure;
hold on;
% set up the gcf and gca before plotting
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'linewidth',1);
set(gca, 'Fontsize', 12);

% Set up title font size
titleFontSize = 28;
axisFontSize = 20;
textFontSize = 16;

% % Add the line at TMS impulse
ylim([-0.02, 0.05]);
y = ylim;
p = plot([0 0], [y(1) y(2)]);
p.Color = 'black';
p.LineWidth = 2;
% Add mark onto the TMS trigger line
t = text(0, y(1)-.005, 'TMS', 'Fontsize', textFontSize);
t.FontWeight = 'bold';
t.HorizontalAlignment = 'center';

% plot the left AMI
p1 = plot(EEG.times, distance_cos_group_ave(1,:), 'lineWidth', 3, 'Color', 'r');
shadedError(EEG.times, distance_cos_group_ave(1,:), distance_cos_group_err(1,:), 'r');
% plot the left UMI
p2 = plot(EEG.times, distance_cos_group_ave(2,:), 'lineWidth', 3, 'Color', 'g');
shadedError(EEG.times, distance_cos_group_ave(2,:), distance_cos_group_err(2,:), 'g');
% plot the right AMI
p3 = plot(EEG.times, distance_cos_group_ave(3,:), 'lineWidth', 3, 'Color', 'b');
shadedError(EEG.times, distance_cos_group_ave(3,:), distance_cos_group_err(3,:), 'b');
% plot the right UMI
p4 = plot(EEG.times, distance_cos_group_ave(4,:), 'lineWidth', 3, 'Color', 'y');
shadedError(EEG.times, distance_cos_group_ave(4,:), distance_cos_group_err(4,:), 'y');
% add the legend
legend([p1, p2, p3, p4], group_name);

ylabel("Cosine Weighted Average of Mahalanobis Distance", 'Fontsize', axisFontSize);
xlabel("Time across the TMS Interval(ms)", 'Fontsize', axisFontSize);
title("Mahalanobis Distances across timepoints on TMS1(spTMS)(Contraposterior Channels)", 'Fontsize', titleFontSize); 
