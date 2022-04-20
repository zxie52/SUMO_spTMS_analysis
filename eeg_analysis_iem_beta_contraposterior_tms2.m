clear;
close all;

%% Start the Journey

% this code using all correct trials instead of supertrials(average 3 trials into 1) on all trials
% for the supertrials usage, please go to the directory: ./old_results_beforeCNS2022

% We have finished the preprocessing

% In this part, I will run the iem on voltage in tms1 through processed data

% This script requires EEGLAB v2020.0 in MATLAB

% The flow for this part:
%                                Loading EEG sets
%                                       |
%                              Prepare for the IEM
%                                       |
%                                 IEM on beta power
%                                       |
%                                 Group Analysis  
%                                       |
%                                    Save out


%% Load the Basics

addpath(genpath('/afs/crc.nd.edu/group/roselab/vol2/zx/new_results/analysis_code'));
addpath(genpath('/afs/crc.nd.edu/group/roselab/vol2/zx/matlab_envi/'));
fpath8 = '/afs/crc.nd.edu/group/roselab/vol2/zx/new_results/eeg_before_IEM';
fpath13 = '/afs/crc.nd.edu/group/roselab/vol2/zx/new_results/IEM_results_beta_tms2';

%%% Warning: due to the warning from iem function("Warning: rank deficient, rank = 6") on SUMO 3001
%%% and SUMO 3017, we have to skip those two subjects first.
%%% Warning: exclude the SUMO 120 because there is no enough corret trials in bin2 for left AMI
%%% switch condition

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_3015' };
% , 'SUMO_3001', 'SUMO_3017'     
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

group_name = {'Left AMI stay', 'Left AMI switch', 'Left UMI stay', 'Left UMI switch', ...
              'Right AMI stay', 'Right AMI switch', 'Right UMI stay', 'Right UMI switch'};

%% Starting IEM for the TMS2 

for l = 1:length(subject)
    for t = 6 % This time we only preprocess the TMS 2 data
        %% Step 1: Load EEG data
        clearvars -except keep fpath8 & fpath13 & subject & l & type & group_name & t;
        
        cd(fpath8);
        EEG = pop_loadset('filename', strcat(subject{l}, '_before_iem_', type{t}, '.set'), 'filepath', fpath8);
        
        EEG_data_doubled = double(EEG.data);% change from single dataset to double
        p = struct2table(EEG.event); % change the EEG.event to a cell array and then categorize into different bins
        
        %% Step 2: Filter trials into four groups: left_AMI, left_UMI, right_AMI, right_UMI
        eventBins = [p.trialnum, p.targetlocation, p.targetlocation2, p.targetorient2, p.epoch, p.leftori, p.rightori];
        % Parameters for later IEM
        nbins = 7;
        binedges = linspace(1,181,nbins+1);
        bincent = round(mean([binedges(1:end-1);binedges(2:end)]));
        
        % filter only correct trials in the probe2
        q = p((p.response2 == 1), :);

        [left_bin_AMI_stay, left_bin_AMI_switch, ...
         left_bin_UMI_stay, left_bin_UMI_switch, ...
         right_bin_AMI_stay, right_bin_AMI_switch, ...
         right_bin_UMI_stay, right_bin_UMI_switch] = filter_bins_for_iem_tms2(q, binedges);
        
        %% IEM
        
        groups = {left_bin_AMI_stay, left_bin_AMI_switch, ...
                  left_bin_UMI_stay, left_bin_UMI_switch, ...
                  right_bin_AMI_stay, right_bin_AMI_switch, ...
                  right_bin_UMI_stay, right_bin_UMI_switch};
        pow_base = [];

        % Run the dothewave function to have pow matrix
        % We can change from [2 50] to [8 13](beta power)or [14 25](beta power)
        % This time we pick the pre-TMS interval as baseline for FFT 
        [pow, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [14 25], 12, 4, 1, [], EEG.times);
  

        %% Extra Step: remove the pre_stim baseline
        % Load the baseline from pre-stimulus period
        cd(fpath8);
        load(strcat('Pre_stimulus_power_baseline_beta', subject{l}, '.mat'))
        tmp = [];

        % Only Pick the beta frequencies(first 6 frequency band)
        tmp = 100 * bsxfun(@rdivide, bsxfun(@minus, pow, pow_base), pow_base);

        % Average the 6 frequencies into one
        beta_power = squeeze(mean(tmp,2));   

        %% Running IEM
        for i = 1 : length(groups) % notice that we used the parallel envi in permutation            
            % have the EEG data and stimlabels for the iem function
            % the output h: epoch# * stimlabels
            h = nonzeros(groups{i});
            [~, colIdcs] = find(groups{i} ~= 0);
            h(:,2) = bincent(colIdcs);
            
            stimlabels = h(:,2);
            super_charge = beta_power(:,:,h(:,1));
            
            % for all posterior electrodes
            % impchan = [23,53,22,54,21,51,19,50,20,48,49,18,10,41,11,42,12,15,44,14,43,45,46,16];
            % for contraposterior channels
            switch i
                case {1 2 3 4}
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
            
            %% Saving and Making the Figures
            
            avgacrosstrials = [];
            avgacrosstrials_perm = [];
            tmp2 = [];
            ChanrespsL = [];
            ChanrespsR = [];
            
            % separate the left and right stimulus
            group_left = 1 : 4;% when i = 1 : 4, we are reconstructing left stimuli
            
            if ismember(i, group_left) % for left stimuli
                avgacrosstrials = mean(chanresp,3); %create a variable that is now 2D
                avgacrosstrials_perm = mean(chanresp_perm,3); %create a variable that is now 2D
                
                % checking if the subject miss the #3 channel in EEG.data
                % for subjects who do not have the 3rd channel, we copy the 5th channel to the 3rd
                if isempty(find(h(:,2) == bincent(3),1))
                    % we made changes on this part from Morgan's original code
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
                
                tmp = mean(chanrespsL,1);
                ChanrespsL = squeeze(tmp);
                
                %plot the encoding model on specific time windows with banded lines
                tw = dsearchn(dstimes',[-200 0]'); %baseline, k = gray
                figure; boundedline(bincent, nanmean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2)))./sqrt(7)),'c','alpha');
                
                tw = dsearchn(dstimes',[20 300]'); %r = red
                hold on; boundedline(bincent, mean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2))))./sqrt(7),'r','alpha');
                
                tw = dsearchn(dstimes',[300 600]'); %m = magenta
                hold on; boundedline(bincent, mean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2))))./sqrt(7) ,'m','alpha');
                
                tw = dsearchn(dstimes',[600 800]'); %b = blue
                hold on; boundedline(bincent, mean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2))))./sqrt(7) ,'b','alpha');
                
                % plot the permutations
                hold on; boundedline(bincent, mean(mean(chanresp_permsL(:,:,34:70),1),3),std(squeeze(chanresp_permsL(:,34:70)))./sqrt(30)*1.96 ,'k','alpha');
                
                title(strcat('Reconstruction from Beta-Power in Contra-Posterior Electrodes', group_name{i}));
                xlabel('Centered Orientation Channel');
                ylabel('Channel Response');
                legend({'error', '-200-0ms','error','20-300ms','error','300-600ms','error','600-800ms'})
                
                cd(fpath13);
                saveas(gcf, strcat('SUMO_', subject{l}, '_iem_', group_name{i}));
                saveas(gcf, strcat('SUMO_', subject{l}, '_iem_', group_name{i}, '.png'));   
                save(strcat('IEM_Exp1_allposterior', subject{l}, '_epoch2_', group_name{i}), 'dstimes', 'stimlabels', 'chanrespsL', 'chanresp_permsL', 'h');
                
            else % for right stimuli
                avgacrosstrials = mean(chanresp,3); %create a variable that is now 2D
                avgacrosstrials_perm = mean(chanresp_perm,3); %create a variable that is now 2D
                
                if isempty(find(h(:,2) == bincent(3),1)) % Should it be 3rd-5th bin shift rather than 1st-6th bin?
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
                tw = dsearchn(dstimes',[-200 0]'); %baseline, k = gray
                figure; boundedline(bincent, nanmean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2)))./sqrt(7)),'c','alpha'); %why is it divided by sqrt(7), not # of trials
                
                tw = dsearchn(dstimes',[20 300]'); %r = red
                hold on; boundedline(bincent, mean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2))))./sqrt(7),'r','alpha');
                
                tw = dsearchn(dstimes',[300 600]'); %m = magenta
                hold on; boundedline(bincent, mean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2))))./sqrt(7) ,'m','alpha');
                
                tw = dsearchn(dstimes',[600 800]'); %b = blue
                hold on; boundedline(bincent, mean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2))))./sqrt(7) ,'b','alpha');
                
                %plot the permutations % why is it 34 : 70, and why is it sqrt(30) * 1.96
                hold on; boundedline(bincent, mean(mean(chanresp_permsR(:,:,34:70),1),3),std(squeeze(chanresp_permsR(:,34:70)))./sqrt(30)*1.96 ,'k','alpha');
                
                title(strcat('Reconstruction from Beta-Power in Contra-Posterior Electrodes', group_name{i}));
                xlabel('Centered Orientation Channel');
                ylabel('Channel Response');
                legend({'error', '-200-0ms','error','20-300ms','error','300-600ms','error','600-800ms'})
                
                cd(fpath13);
                saveas(gcf, strcat('SUMO_', subject{l}, '_iem_', group_name{i}));
                saveas(gcf, strcat('SUMO_', subject{l}, '_iem_', group_name{i}, '.png'));   
                save(strcat('IEM_Exp1_allposterior', subject{l}, '_epoch2_', group_name{i}), 'dstimes', 'stimlabels', 'chanrespsR', 'chanresp_permsR', 'h');
            end
        end
    end
end

%% Merging the IEM for group analysis(both left and right AMI)
clearvars -except keep fpath8 & fpath13 & subject & l & group_name;
cd(fpath13);

nbins = 7;
binedges = linspace(1,181,nbins+1);
bincent = round(mean([binedges(1:end-1);binedges(2:end)]));

% Load Parameters
centerind = 4;
nchan = 7;

% Set up title font size
titleFontSize = 25;
axisFontSize = 22;
        

for i = 1 : length(group_name)
    group_left = 1 : 4;% when i = 1 : 4, we are reconstructing left stimuli
    temp1_CR1L = [];
    temp1_CRP1L = [];
    temp1_CR1R = [];
    temp1_CRP1R = [];
    
    for j = 1: length(subject)
        load(strcat('IEM_Exp1_allposterior', subject{j}, '_epoch2_', group_name{i}));
     
        if ismember(i, group_left) % for left stimuli
            %squeeze the channel response
            tmp1 = mean(chanrespsL,1);
            chanrespsL = squeeze(tmp1);
            
            tmp3 = mean(chanresp_permsL,1);
            chanresp_permsL = squeeze(tmp3);
            
            temp1_CR1L(j,:,:) = chanrespsL;
            temp1_CRP1L(j,:,:) = chanresp_permsL;
            
        else % for right stimuli
            tmp2 = mean(chanrespsR,1);
            chanrespsR = squeeze(tmp2);
            
            tmp4 = mean(chanresp_permsR,1);
            chanresp_permsR = squeeze(tmp4);
            
            temp1_CR1R(j,:,:) = chanrespsR;
            temp1_CRP1R(j,:,:) = chanresp_permsR;
        end
    end % for nine subjects
    
    chanrespsL_allsub = temp1_CR1L;
    chanrespsR_allsub = temp1_CR1R;
    chanresp_permsL_allsub = temp1_CRP1L;
    chanresp_permsR_allsub = temp1_CRP1R;
    
    save(strcat('IEM_Exp1_contraposterior_allsubjects_epoch2_', group_name{i}, '_iem'), ...
    'dstimes', 'stimlabels', 'chanrespsL_allsub', 'chanresp_permsL_allsub', ...
    'chanrespsR_allsub', 'chanresp_permsR_allsub'); 

%% Plotting Left 
    % plotting the group level is different from plotting the single subject
    if ismember(i, group_left) % for left stimuli
        for a = 1 : size(chanrespsL_allsub, 1) % # of subjects
            for b = 1 : size(chanrespsL_allsub, 2) % # of bins
                % The channel response is 1000ms (250 time points)
                % Use 12 as the time window to get 50ms per sample
                condL(a,b,:) = smoothdata(chanrespsL_allsub(a,b,:), 'gaussian', 12);
            end
        end
        
        ChanrespsL_allsub = squeeze(nanmean(condL,1));
        
        % plot the encoding model on specific time windows with banded lines(left)
        % Generally, boundedline(bincent, mean of the timepoint on bins(7 * 1 double), ...
        %                           std on the subjects(1st dimension), of the mean on the timepoint(3rd dimension))
        tw = dsearchn(dstimes',[-200 0]'); %baseline, k = gray
        figure; boundedline(bincent, squeeze(nanmean(ChanrespsL_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsL_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'c','alpha');

        tw = dsearchn(dstimes',[20 300]'); %r = red
        boundedline(bincent, squeeze(nanmean(ChanrespsL_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsL_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'r','alpha');

        tw = dsearchn(dstimes',[300 600]'); %m = magenta
        boundedline(bincent, squeeze(nanmean(ChanrespsL_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsL_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'m','alpha');

        tw = dsearchn(dstimes',[600 800]'); %b = blue
        boundedline(bincent, squeeze(nanmean(ChanrespsL_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsL_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'b','alpha');

        % plot the permutations % why is it 34 : 70, and why is it sqrt(30) * 1.96
        % Generally, boundedline(bincent, mean of the subjects of the mean of timepoint on bins(7 * 1 double), ...
        %                           std on the subjects(1st dimension), of the mean on the timepoint(3rd dimension))
        boundedline(bincent, squeeze(mean(nanmean(chanresp_permsL_allsub(:,:,34:70),3),1)),squeeze(std(mean(chanresp_permsL_allsub(:,:,34:70),3)))./sqrt(30)*1.96,'k','alpha');

        hold on;
        ylim([-.2 .6]);

        title(strcat("Reconstruction from Beta-Power in Contra-Posterior Electrodes (TMS2 ", group_name{i}, ")"));
        xlabel("Centered Orientation Channel");
        ylabel("Channel Response");
        legend({'error', '-200-0ms','error','20-300ms','error','300-600ms','error','600-800ms'})

        saveas(gcf, strcat('All_Subjects_iem_', group_name{i},'.png') );    

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
        for k = 1: length(S)
            if S(k).Area >= threshSize
                idx{j} = S(k).PixelIdxList';
                j = j + 1;
            end
        end
        % return the idx for multiple clusters
        
        %% Plotting heatmap
        % the process for generating the heat map
        orioffset = linspace(-90, 90, nchan); %for plotting
        figure; contourf(dstimes, orioffset, squeeze(mean(condL, 1)),30,'linec','none'); %make a contour plot
        xlim([min(dstimes) max(dstimes)]); %can trim the time range of the plot here if you want. Right now setting to min and max does nothing.
        h = colorbar;
        caxis([-0.1 0.5]); %scale the colorbar axis to just region of interest (-0.1 to 0.6)

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
        
        ylabel(h, 'Channel Response (µV^2)','Fontsize', axisFontSize)
        title(strcat("Reconstruction from Beta-Power in Contra-Posterior Electrodes (TMS2 ", group_name{i}, ")"), 'Fontsize', titleFontSize)
        xlabel('Time from TMS onset(msec)','Fontsize', axisFontSize)
        ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)
        savefig(strcat('All_Subjects_iem_', group_name{i},'_heatmap'));
        saveas(gcf, strcat('All_Subjects_iem_', group_name{i},'_heatmap.png'));
    else
        %% Plotting Right
        for a = 1 : size(chanrespsR_allsub, 1) % # of subjects
            for b = 1 : size(chanrespsR_allsub, 2) % # of bins
                % The channel response is 1000ms (250 time points)
                % Use 12 as the time window to get 50ms per sample
                condR(a,b,:) = smoothdata(chanrespsR_allsub(a,b,:), 'gaussian', 12);
            end
        end
        
        ChanrespsR_allsub = squeeze(nanmean(condR,1));
        %plot the encoding model on specific time windows with banded lines(right)
        tw = dsearchn(dstimes',[-200 0]'); %baseline, k = gray
        figure; boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'c','alpha');

        tw = dsearchn(dstimes',[20 300]'); %r = red
        boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'r','alpha');

        tw = dsearchn(dstimes',[300 600]'); %m = magenta
        boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'m','alpha');

        tw = dsearchn(dstimes',[600 800]'); %b = blue
        boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'b','alpha');

        % plot the permutations % why is it 34 : 70, and why is it sqrt(30) * 1.96
        boundedline(bincent, squeeze(mean(nanmean(chanresp_permsR_allsub(:,:,34:70),3),1)),squeeze(std(mean(chanresp_permsR_allsub(:,:,34:70),3)))./sqrt(30)*1.96,'k','alpha');
        
        hold on;
        ylim([-.2 .6]);

        title(strcat('Reconstruction from Beta-Power in Contra-Posterior Electrodes(TMS2 ', group_name{i}, ')'));
        xlabel('Centered Orientation Channel');
        ylabel('Channel Response');
        legend({'error', '-200-0ms','error','20-300ms','error','300-600ms','error','600-800ms'})
        
        saveas(gcf, strcat('All_Subjects_iem_', group_name{i}, '.png')); 
        
        
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
        for k = 1: length(S)
            if S(k).Area >= threshSize
                idx{j} = S(k).PixelIdxList';
                j = j + 1;
            end
        end
        
        %% Plotting heatmap
        % the process for generating the heat map
        orioffset = linspace(-90, 90, nchan); %for plotting
        figure; contourf(dstimes, orioffset, squeeze(mean(condR, 1)),30,'linec','none'); %make a contour plot
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
        
        ylabel(h, 'Channel Response (µV^2)','Fontsize', axisFontSize)
        title(strcat("Reconstruction from Beta-Power in Contra-Posterior Electrodes (TMS2 ", group_name{i}, ")"), 'Fontsize', titleFontSize)
        xlabel('Time from TMS onset(msec)','Fontsize', axisFontSize)
        ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)
        savefig(strcat('All_Subjects_iem_', group_name{i},'_heatmap'));
        saveas(gcf, strcat('All_Subjects_iem_', group_name{i},'_heatmap.png'));
    end
end

% finished
load handel;
sound(y, Fs);