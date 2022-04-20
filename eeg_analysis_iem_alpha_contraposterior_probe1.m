close all;
clear;

%% Start the Journey

% this code using all correct trials instead of supertrials(average 3 trials into 1) on all trials
% for the supertrials usage, please go to the directory: ./old_results_beforeCNS2022

% This script requires EEGLAB v2020.0 in MATLAB

% The flow for this part:
%                                Loading EEG sets
%                                       |
%                              Prepare for the IEM
%                                       |
%                                 IEM on voltage
%                                       |
%                                 Group Analysis  
%                                       |
%                                    Save out

%% Load the Basics

addpath(genpath('E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\analysis_code'));
fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120',...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};
        
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

group_name = {'left probe', 'right probe'};

%% Picking the IEM on target orientation or test orientation
condition = input("Do you want the IEM on testorient(test orientation) or targetorient(target orientation)?\n", 's');
while ~ismember(condition, {'testorient', 'targetorient'})
    warning("The input condition is not available, please re-type the condition.");
    condition = input("Do you want the IEM on testorient(test orientation) or targetorient(target orientation)?\n", 's');
end

if strcmp(condition, 'testorient')
    fpath13 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\IEM_results_alpha_probe1_testorient';
end

if strcmp(condition, 'targetorient')
    fpath13 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\IEM_results_alpha_probe1_targetorient';
end

%% Start the IEM
for l = 1:length(subject)
    %1 : length(subject)
    %for t = 1 : length(type) % Picking different types
    for t = 4 % This time we only preprocess the Probe1 data
        %% Step 1: Load EEG data
        clearvars -except keep fpath8 & fpath13 & subject & l & type & group_name & t & condition;
        
        cd(fpath8);
        EEG = pop_loadset('filename', strcat(subject{l}, '_before_iem_', type{t}, '.set'), 'filepath', fpath8);
        
        EEG_data_doubled = double(EEG.data);% change from single dataset to double
        p = struct2table(EEG.event); % change the EEG.event to a cell array and then categorize into different bins 
       
        if strcmp(condition, 'testorient')
            % from the test script, p.accuracy = p.testorient - p.targetorient 
            for i = 1 : height(p)
                p.testorient(i) = p.targetorient(i) + p.accuracy(i);
            end

            % for probe orientations which are negative, turn it to 180 + testorient
            % for probe orientations which are above 180. adjust it to testorient - 180
            for i = 1 : height(p)
                if p.testorient(i) > 180
                    p.testorient(i) = p.testorient(i) - 180;
                elseif p.testorient(i) < 0
                        p.testorient(i) = 180 + p.testorient(i);
                end
            end
            % move the testorient column right after the targetorient column
            p = movevars(p, 'testorient', 'After', 'targetorient');
        end 
        
        %% Step 3: filter trials into two groups: left_probe, right_probe
        %get the event trial number; target location1; test orientation; epoch#; left orientation(stimulus); right orientaion(stimulus); 
        
        %Filter trials into left and right group
        %bins(1.trial#, 2.trialloc, 3.targetori, 4.epoch#, 5.leftori, 6.rightori)
        
        nbins = 7;
        binedges = round(linspace(1,181,nbins+1));
        bincent = round(mean([binedges(1:end-1);binedges(2:end)]));
        
        % filter only correct trials in the probe1
        q = p((p.response == 1), :);

        [leftTrial, leftStimlabels, rightTrial, rightStimlabels] = filter_bins_for_iem_probe(condition, type{t}, p, binedges, bincent);
        
        %% Step 4: IEM
        groups = {leftTrial, rightTrial};
        
        %% Extra Step: run the FFT on EEG data and subtract the pre-stimulus baseline
        pow_base = [];

        % Run the dothewave function to have pow matrix
        % We can change from [2 50] to [8 13](alpha power)
        % This time we pick the pre-stimulus interval as baseline for FFT 
        [pow, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [8 13], 6, 4, 1, [], EEG.times);

        % Load the baseline from pre-stimulus period
        cd(fpath8);
        load(strcat('Pre_stimulus_power_baseline', subject{l}, '.mat'))
        tmp = [];
        
        % Only Pick the alpha frequencies(first 6 frequency band)
        tmp = 100 * bsxfun(@rdivide, bsxfun(@minus, pow, pow_base), pow_base);

        % Average the 6 frequencies into one
        alpha_power = squeeze(mean(tmp,2));   

        %% Running the IEM on the probe
        for i = 1 : length(group_name) % notice that we used the parallel envi in permutation
            chanrespsL = [];
            chanrespsR = [];
            chanresp_permsL = [];
            chanresp_permsR = [];
            
            % have the EEG data and stimlabels for the iem function
            % the output h: epoch# * stimlabels
            h = nonzeros(groups{i});
            [~, colIdcs] = find(groups{i} ~= 0);
            h(:,2) = bincent(colIdcs);
            
            stimlabels = h(:,2);
            super_charge = alpha_power(:,:,h(:,1));

            nbins = 7;
            binedges = linspace(1,181,nbins+1);
            bincent = round(mean([binedges(1:end-1);binedges(2:end)]));

            clear chanresp & weights & dstimes & chanresp_perm;

            % for all posterior electrodes
            % impchan = [23,53,22,54,21,51,19,50,20,48,49,18,10,41,11,42,12,15,44,14,43,45,46,16]; 

            %for contraposterior electrodes
            switch i
                case 1 % for left probes
                    ROI = {'CP2', 'CP4', 'CP6', 'TP8', 'TP10', ...
                           'P2', 'P4', 'P6', 'P8', ...
                           'PO4', 'PO8', 'O2'};
                    impchan = find(ismember({EEG.chanlocs.labels}, ROI)); %channels in R hem
                otherwise % for right probes
                    ROI = {'CP1', 'CP3', 'CP5', 'TP7', 'TP9', ...
                           'P1', 'P3', 'P5', 'P7', ...
                           'PO3', 'PO7', 'O1'};
                    impchan = find(ismember({EEG.chanlocs.labels}, ROI)); %channels in L hem
            end

            % IEM part
            [chanresp, ~, dstimes] = iemori(super_charge(impchan,:,:),stimlabels,4,EEG.times);

            % permutation test to get null hypothesis - takes a LONG time
            % in the permutation, the stimlabels need to be shuffled
            nperm = 100;
            parfor p = 1:nperm
                disp(p)
                [tmp, ~, ~] = iemori(super_charge(impchan,:,:),stimlabels(randperm(length(stimlabels))),4,EEG.times); 
                %the 4,EEG.times is the vector of times with a downsampling factor of 4
                chanresp_perm(:,:,p) = mean(tmp,3);
            end
            
            %% Step 5: Saving and Making Individual Figures
            avgacrosstrials = [];
            avgacrosstrials_perm = [];
            tmp2 = [];
            ChanrespsL = [];
            ChanrespsR = [];

            % Separate the left and right circumstances
            if i == 1
                avgacrosstrials = mean(chanresp,3); %create a variable that is now 2D
                avgacrosstrials_perm = mean(chanresp_perm,3); %create a variable that is now 2D
            
                % checking if the subject miss the #3 channel in EEG.data
                % for subjects who do not have the 3rd channel, we copy the 5th channel to the 3rd
                if isempty(find(h(:,2) == bincent(3),1)) 
                    chanrespsL(1,3,:) = avgacrosstrials(4,:);
                    chanrespsL(1,1:2,:) = avgacrosstrials(1:2,:);
                    chanrespsL(1,4:7,:) = avgacrosstrials(3:6,:);
                    chanresp_permsL(1,3,:) = avgacrosstrials_perm(4,:);
                    chanresp_permsL(1,1:2,:) = avgacrosstrials_perm(1:2,:);
                    chanresp_permsL(1,4:7,:) = avgacrosstrials_perm(3:6,:);
                else
                    chanrespsL(1,:,:) = mean(chanresp,3);
                    chanresp_permsL(1,:,:) = mean(chanresp_perm,3);
                end

                tmp2 = mean(chanrespsL,1);
                ChanrespsL = squeeze(tmp2);
                
                %plot the encoding model on specific time windows with banded lines
                tw = dsearchn(dstimes',[-100 0]'); %baseline, k = gray
                figure; boundedline(bincent, nanmean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2)))./sqrt(7)),'c','alpha'); 

                tw = dsearchn(dstimes',[20 200]'); %r = red
                hold on; boundedline(bincent, mean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2))))./sqrt(7),'r','alpha');

                tw = dsearchn(dstimes',[200 400]'); %m = magenta
                hold on; boundedline(bincent, mean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2))))./sqrt(7) ,'m','alpha');
                
                % plot the permutations 
                hold on; boundedline(bincent, mean(mean(chanresp_permsL(:,:,34:70),1),3),std(squeeze(chanresp_permsL(:,34:70)))./sqrt(30)*1.96 ,'k','alpha');
                
                hold on;
                ylim([-.2 .6]);
                
                title(strcat('Reconstruction from Alpha Power in Contra-Posterior Electrodes(', group_name{i}, ')'));
                xlabel('Centered Orientation Channel');
                ylabel('Channel Response');
                legend({'error', '-100-0ms','error','20-200ms','error','200-400ms'})

                cd(fpath13);
                saveas(gcf, strcat('SUMO_', subject{l}, '_iem_', group_name{i}, '.png'));    
                save(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{l}, '_epoch1_', group_name{i}), 'dstimes', 'stimlabels', 'chanrespsL', 'chanresp_permsL', 'h'); 
            else % for right stimuli
                avgacrosstrials = mean(chanresp,3); %create a variable that is now 2D
                avgacrosstrials_perm = mean(chanresp_perm,3); %create a variable that is now 2D
            
                % checking if the subject miss the #3 channel in EEG.data
                % for subjects who do not have the 3rd channel, we copy the 5th channel to the 3rd
                if isempty(find(h(:,2) == bincent(3),1)) 
                    chanrespsR(1,3,:) = avgacrosstrials(4,:);
                    chanrespsR(1,1:2,:) = avgacrosstrials(1:2,:);
                    chanrespsR(1,4:7,:) = avgacrosstrials(3:6,:);
                    chanresp_permsR(1,3,:) = avgacrosstrials_perm(4,:);
                    chanresp_permsR(1,1:2,:) = avgacrosstrials_perm(1:2,:);
                    chanresp_permsR(1,4:7,:) = avgacrosstrials_perm(3:6,:);
                else
                    chanrespsR(1,:,:) = mean(chanresp,3);
                    chanresp_permsR(1,:,:) = mean(chanresp_perm,3);
                end

                tmp2 = mean(chanrespsR,1);
                ChanrespsR = squeeze(tmp2);
                
                %plot the encoding model on specific time windows with banded lines
                tw = dsearchn(dstimes',[-100 0]'); %baseline, k = gray
                figure; boundedline(bincent, nanmean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2)))./sqrt(7)),'c','alpha'); 

                tw = dsearchn(dstimes',[20 200]'); %r = red
                hold on; boundedline(bincent, mean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2))))./sqrt(7),'r','alpha');

                tw = dsearchn(dstimes',[200 400]'); %m = magenta
                hold on; boundedline(bincent, mean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2))))./sqrt(7) ,'m','alpha');

                % plot the permutations 
                hold on; boundedline(bincent, mean(mean(chanresp_permsR(:,:,34:70),1),3),std(squeeze(chanresp_permsR(:,34:70)))./sqrt(30)*1.96 ,'k','alpha');
                
                hold on;
                ylim([-.2 .6]);
                
                title(strcat('Reconstruction from Alpha Power in Contra-Posterior Electrodes(', group_name{i}, ')'));
                xlabel('Centered Orientation Channel');
                ylabel('Channel Response');
                legend({'error', '-100-0ms','error','20-200ms','error','200-400ms'})

                cd(fpath13);
                saveas(gcf, strcat('SUMO_', subject{l}, '_iem_', group_name{i}, '.png'));    
                save(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{l}, '_epoch1_', group_name{i}), 'dstimes', 'stimlabels', 'chanrespsR', 'chanresp_permsR', 'h'); 
            end
        end
    end
end

%% Merging the IEM for group analysis(both left and right probe)
clearvars -except keep fpath8 & fpath13 & subject & group_name & condition;
cd(fpath13);

for i = 1:length(subject)
    for j = 1:2
        load(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{i}, '_epoch1_', group_name{j}, '.mat'));
    end
    
%     %squeeze the channel response
%     chanrespsL = squeeze(chanrespsL);
%     chanrespsR = squeeze(chanrespsR);
%     chanresp_permsL = squeeze(chanresp_permsL);
%     chanresp_permsR = squeeze(chanresp_permsR);

    temp1_CR1L(i,:,:) = chanrespsL;
    temp1_CR1R(i,:,:) = chanrespsR;
    temp1_CRP1L(i,:,:) = chanresp_permsL;
    temp1_CRP1R(i,:,:) = chanresp_permsR;

end

chanrespsL_allsub = temp1_CR1L;
chanrespsR_allsub = temp1_CR1R;
chanresp_permsL_allsub = temp1_CRP1L;
chanresp_permsR_allsub = temp1_CRP1R;

save(strcat('IEM_Exp1_contraposterior_allsubjects_epoch1_probe_iem.mat'), 'dstimes', 'stimlabels', 'chanrespsL_allsub', 'chanresp_permsL_allsub', 'chanrespsR_allsub', 'chanresp_permsR_allsub'); 

%% Plotting Left 
% plotting the group level is different from plotting the single subject
nbins = 7;
binedges = linspace(1,181,nbins+1);
bincent = round(mean([binedges(1:end-1);binedges(2:end)]));

ChanrespsL_allsub = squeeze(nanmean(chanrespsL_allsub,1));

% plot the encoding model on specific time windows with banded lines(left)
% Generally, boundedline(bincent, mean of the timepoint on bins(7 * 1 double), ...
%                           std on the subjects(1st dimension), of the mean on the timepoint(3rd dimension))
tw = dsearchn(dstimes',[-100 0]'); %baseline, k = gray
figure; boundedline(bincent, squeeze(nanmean(ChanrespsL_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsL_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'c','alpha');

tw = dsearchn(dstimes',[20 200]'); %r = red
boundedline(bincent, squeeze(nanmean(ChanrespsL_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsL_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'r','alpha');

tw = dsearchn(dstimes',[200 400]'); %m = magenta
boundedline(bincent, squeeze(nanmean(ChanrespsL_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsL_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'m','alpha');

% plot the permutations % why is it 34 : 70, and why is it sqrt(30) * 1.96
% Generally, boundedline(bincent, mean of the subjects of the mean of timepoint on bins(7 * 1 double), ...
%                           std on the subjects(1st dimension), of the mean on the timepoint(3rd dimension))

boundedline(bincent, squeeze(mean(nanmean(chanresp_permsL_allsub(:,:,34:70),3),1)),squeeze(std(mean(chanresp_permsL_allsub(:,:,34:70),3)))./sqrt(30)*1.96,'k','alpha');

hold on;
ylim([-.2 .6]);
if strcmp(condition, 'testorient')
    title('IEM Left Stimulus Reconstruction from Contraposterior Electrodes(Probe1 testorient)');
else
    title('IEM Left Stimulus Reconstruction from Contraposterior Electrodes(Probe1 targetorient)');
end
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')
legend({'error', '-100-0ms','error','20-200ms','error','200-400ms'})

xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')

saveas(gcf, 'allsubjects_iem_left_probe.png');
%% Plotting Right

ChanrespsR_allsub = squeeze(nanmean(chanrespsR_allsub,1));
%plot the encoding model on specific time windows with banded lines(right)
tw = dsearchn(dstimes',[-100 0]'); %baseline, k = gray
figure; boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'c','alpha');

tw = dsearchn(dstimes',[20 200]'); %r = red
boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'r','alpha');

tw = dsearchn(dstimes',[200 400]'); %m = magenta
boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'m','alpha');

% plot the permutations % why is it 34 : 70, and why is it sqrt(30) * 1.96
boundedline(bincent, squeeze(mean(nanmean(chanresp_permsR_allsub(:,:,34:70),3),1)),squeeze(std(mean(chanresp_permsR_allsub(:,:,34:70),3)))./sqrt(30)*1.96,'k','alpha');

hold on;
ylim([-.2 .6]);

if strcmp(condition, 'testorient')
    title('IEM Right Stimulus Reconstruction from Contraposterior Electrodes(Probe1 testorient)');
else
    title('IEM Right Stimulus Reconstruction from Contraposterior Electrodes(Probe1 targetorient)');
end
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')
legend({'error', '-100-0ms','error','20-200ms','error','200-400ms'})

xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')

saveas(gcf, 'allsubjects_iem_right_probe.png');

%% Plottingg the heatmaps for the left and right probe
% Load Parameters
centerind = 4;
nchan = 7;

% Set up title font size
titleFontSize = 25;
axisFontSize = 22;

% Smooth the data to 50ms window by Gaussian filter
for k =  1 : 2
    if k == 1
        for i = 1 : size(chanrespsL_allsub, 1) % # of subjects
            for j = 1 : size(chanrespsL_allsub, 2) % # of bins
                % The channel response is 1000ms (250 time points)
                % Use 12 as the time window to get 50ms per sample
                chanrespsL_allsub(i,j,:) = smoothdata(chanrespsL_allsub(i,j,:), 'gaussian', 12);
            end
        end
        condL = chanrespsL_allsub;
        
        %% Do the linear fitting function for IEM channel response and permutation
        [slope1, slope2] = linFitIEM(condL, chanresp_permsL_allsub);

        %% Do the t-test between the channel response on the center channel and the permutation
        % do the t-test, right tailed, .05
        [~, p] = ttest(slope1,slope2,.05,'right');

        %% Do the Cluster-based Permutation to correct the p-value
        % Use cluster based permutation to corrrect the p-value
        threshSize = clustthresh1D(slope1, slope2, 10000, 'right');
        S=regionprops(p<.05,'PIxelIdxList','Area');
        idx = {};
        j = 1;
        for i = 1: length(S)
            if S(i).Area >= threshSize
                idx{j} = S(i).PixelIdxList';
                j = j + 1;
            end
        end
        % return the idx for multiple clusters
    
    else % for the right stimulus
        for i = 1 : size(chanrespsR_allsub, 1) % # of subjects
            for j = 1 : size(chanrespsR_allsub, 2) % # of bins
                % The channel response is 1000ms (250 time points)
                % Use 12 as the time window to get 50ms per sample
                chanrespsR_allsub(i,j,:) = smoothdata(chanrespsR_allsub(i,j,:), 'gaussian', 12);
            end
        end
        condR = chanrespsR_allsub;
        
        %% Do the linear fitting function for IEM channel response and permutation
        [slope1, slope2] = linFitIEM(condR, chanresp_permsR_allsub);

        %% Do the t-test between the channel response on the center channel and the permutation
        % do the t-test, right tailed, .05
        [~, p] = ttest(slope1,slope2,.05,'right');

        %% Do the Cluster-based Permutation to correct the p-value
        % Use cluster based permutation to corrrect the p-value
        threshSize = clustthresh1D(slope1, slope2, 10000, 'right');
        S=regionprops(p<.05,'PIxelIdxList','Area');
        idx = {};
        j = 1;
        for i = 1: length(S)
            if S(i).Area >= threshSize
                idx{j} = S(i).PixelIdxList';
                j = j + 1;
            end
        end
        % return the idx for multiple clusters  
    end
    
    % the process for generating the heat map
    orioffset = linspace(-90, 90, nchan); %for plotting
    if k == 1
        figure;contourf(dstimes, orioffset, squeeze(mean(condL, 1)),30,'linec','none'); %make a contour plot
    else
        figure;contourf(dstimes, orioffset, squeeze(mean(condR, 1)),30,'linec','none'); %make a contour plot
    end
    xlim([min(dstimes) max(dstimes)]); %can trim the time range of the plot here if you want. Right now setting to min and max does nothing.
    h = colorbar;
    caxis([-0.1 0.5]); %scale the colorbar axis to just region of interest (-0.1 to 0.6)

    % Plotting the uncorrected p-value on the heatmap
    hold on;
    if ~isempty(p(p<.05))
        plot(dstimes(p<.05), 80, 'g.', 'MarkerSize', 25);
    end

    % Plotting the clusters on the heatmap
    hold on;
    if ~isempty(cell2mat(idx))
        plot(dstimes(cell2mat(idx)), 75, 'r.', 'MarkerSize', 25);
    else
        warning("There is no significant p-value after clusted-based permutation correction");
        warning("The threshold size is, %d", threshSize);
        warning("The largest cluster size is, %d", max([S.Area]));
    end

    % Opening the figure in the full size
    set(gcf, 'Position', get(0, 'Screensize'));
    % determine if it is left or right channel response
    if k == 1
        ylabel(h, 'Channel Response (µV^2)','Fontsize', axisFontSize)
        if strcmp(condition, 'testorient')
            title('Reconstruction from Alpha Power in ContraPosterior Channels Left probe(Probe1)(testorient)', 'Fontsize', titleFontSize)
        elseif strcmp(condition, 'targetorient')
            title('Reconstruction from Alpha Power in ContraPosterior Channels Left probe(Probe1)(targetorient)', 'Fontsize', titleFontSize)
        end            
        xlabel('Time from probe onset(msec)','Fontsize', axisFontSize)
        ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)
        savefig('allsubjects_iem_left_probe_heatmap');
        saveas(gcf, 'allsubjects_iem_left_probe_heatmap.png');
    elseif k == 2 
        ylabel(h, 'Channel Response (µV^2)','Fontsize', axisFontSize)
        if strcmp(condition, 'testorient')
            title('Reconstruction from Alpha Power in ContraPosterior Channels Right probe(Probe1)(testorient)', 'Fontsize', titleFontSize)
        elseif strcmp(condition, 'targetorient')
            title('Reconstruction from Alpha Power in ContraPosterior Channels Right probe(Probe1)(targetorient)', 'Fontsize', titleFontSize)
        end   
        xlabel('Time from probe onset(msec)','Fontsize', axisFontSize)
        ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)
        savefig('allsubjects_iem_right_probe_heatmap');
        saveas(gcf, 'allsubjects_iem_right_probe_heatmap.png');
    end
end

% finished
load handel;
sound(y, Fs);