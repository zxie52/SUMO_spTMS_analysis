close all;
clear;

%% Start the Journey

% this code using all correct trials instead of supertrials(average 3 trials into 1) on all trials
% for the supertrials usage, please go to the directory: ./old_results_beforeCNS2022

% We have finished the preprocessing

% In this part, I will run the iem on voltage in cue1 through processed data

% This script requires EEGLAB v2020.0 in MATLAB

% The flow for this part:
%                                Loading EEG sets
%                                       |
%                              Prepare for the IEM
%                                       |
%                                 IEM on alpha-power
%                                       |
%                                 Group Analysis  
%                                       |
%                                    Save out
%% Load the Basics

addpath(genpath("/afs/crc.nd.edu/group/roselab/vol2/zx/new_results/analysis_code"));
addpath(genpath('/afs/crc.nd.edu/group/roselab/vol2/zx/matlab_envi/'));
fpath8 = '/afs/crc.nd.edu/group/roselab/vol2/zx/new_results/eeg_before_IEM';
fpath13 = '/afs/crc.nd.edu/group/roselab/vol2/zx/new_results/IEM_results_alpha_tms1';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120', ...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};
       
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

group_name = {'left AMI', 'left UMI', 'right AMI', 'right UMI'};

%% Start the IEM
for l = 1:length(subject)
    %1 : length(subject)
    %for t = 1 : length(type) % Picking different types
    for t = 3 % This time we only preprocess the TMS 1 data
        %% Step 1: Load EEG data
        clearvars -except keep fpath8 & fpath13 & subject & l & type & group_name & t;
        
        cd(fpath8);
        EEG = pop_loadset('filename', strcat(subject{l}, '_before_iem_', type{t}, '.set'), 'filepath', fpath8);
        
        %Remove the eventmarkers of 'S 50'
        EEG.event(strcmp({EEG.event.type}, 'S 50')) = [];
            
        p = struct2table(EEG.event); % change the EEG.event to a cell array and then categorize into different bins
        
        %% Step 2: filter trials into four groups: left_AMI, left_UMI, right_AMI, right_UMI
        %get the event trial number; trial location1; target orientation; epoch#; left orientation(stimulus); right orientaion(stimulus); 
        eventBins = [p.trialnum, p.targetlocation, p.targetorient, p.epoch, p.leftori, p.rightori];
        
        % we will take 7 bins for the IEM
        nbins = 7;
        % calculate the edge between each bin
        binedges = linspace(1,181,nbins+1);
        % calculate the center of each bin
        bincent = round(mean([binedges(1:end-1);binedges(2:end)]));
        
        % filter only correct trials in the probe1
        q = p((p.response == 1), :);
        
        %Filter trials into left and right group
        [left_bin_AMI, left_bin_UMI, right_bin_AMI, right_bin_UMI] = filter_bins_for_iem_tms(q, binedges, bincent);
        
        %% Step 3: IEM
        groups = {left_bin_AMI, left_bin_UMI, right_bin_AMI, right_bin_UMI};
        pow_base = [];

        % Run the dothewave function to have pow matrix
        % We can change from [2 50] to [8 13](alpha power)
        % This time we pick the pre-TMS interval as baseline for FFT 
        [pow, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [8 13], 6, 4, 1, [], EEG.times);
  

        %% Extra Step: remove the pre_stim baseline
        % Load the baseline from pre-stimulus period
        cd(fpath8);
        load(strcat('Pre_stimulus_power_baseline', subject{l}, '.mat'))
        tmp = [];

        % Only Pick the alpha frequencies(first 6 frequency band)
        tmp = 100 * bsxfun(@rdivide, bsxfun(@minus, pow, pow_base), pow_base);

        % Average the 6 frequencies into one
        alpha_power = squeeze(mean(tmp,2));   
        
        for i = 1 : length(group_name) % notice that we used the parallel envi in permutation
            chanrespsL = [];
            chanrespsR = [];
            chanresp_permsL = [];
            chanresp_permsR = [];
            
            %% Trials for four different stimuli (AMI + UMI)            
            % have the EEG data and stimlabels for the iem function
            % the output h: epoch# * stimlabels
            h = nonzeros(groups{i});
            [~, colIdcs] = find(groups{i} ~= 0);
            h(:,2) = bincent(colIdcs);
            
            stimlabels = h(:,2);
            super_charge = alpha_power(:,:,h(:,1));

            clear chanresp & weights & dstimes & chanresp_perm;

            % for all posterior electrodes
            % impchan = [23,53,22,54,21,51,19,50,20,48,49,18,10,41,11,42,12,15,44,14,43,45,46,16]; 

            %for contraposterior electrodes
            switch i
                case 1
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

            %% IEM part
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
            
            %% Step 4: Saving and Making Individual Figures
            avgacrosstrials = [];
            avgacrosstrials_perm = [];
            tmp2 = [];
            ChanrespsL = [];
            ChanrespsR = [];

            % Separate the left and right circumstances
            if ismember(i, [1 2])
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
                tw = dsearchn(dstimes',[-200 -20]'); %baseline, k = gray
                figure; boundedline(bincent, nanmean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2)))./sqrt(7)),'c','alpha'); 

                tw = dsearchn(dstimes',[30 300]'); %r = red
                hold on; boundedline(bincent, mean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2))))./sqrt(7),'r','alpha');

                tw = dsearchn(dstimes',[300 600]'); %m = magenta
                hold on; boundedline(bincent, mean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2))))./sqrt(7) ,'m','alpha');

                tw = dsearchn(dstimes',[600 800]'); %b = blue
                hold on; boundedline(bincent, mean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2))))./sqrt(7) ,'b','alpha');

                % plot the permutations 
                hold on; boundedline(bincent, mean(mean(chanresp_permsL(:,:,34:70),1),3),std(squeeze(chanresp_permsL(:,34:70)))./sqrt(30)*1.96 ,'k','alpha');

                title(strcat('Reconstruction from Alpha Power in Contra-Posterior Electrodes(', group_name{i}, ')'));
                xlabel('Centered Orientation Channel');
                ylabel('Channel Response');
                legend({'error', '-200-0ms','error','30-300ms','error','300-600ms','error','600-800ms'})

                cd(fpath13);
                saveas(gcf, strcat('SUMO_', subject{l}, '_iem_', group_name{i}));
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
                tw = dsearchn(dstimes',[-200 -20]'); %baseline, k = gray
                figure; boundedline(bincent, nanmean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2)))./sqrt(7)),'c','alpha'); 

                tw = dsearchn(dstimes',[30 300]'); %r = red
                hold on; boundedline(bincent, mean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2))))./sqrt(7),'r','alpha');

                tw = dsearchn(dstimes',[300 600]'); %m = magenta
                hold on; boundedline(bincent, mean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2))))./sqrt(7) ,'m','alpha');

                tw = dsearchn(dstimes',[600 800]'); %b = blue
                hold on; boundedline(bincent, mean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2))))./sqrt(7) ,'b','alpha');

                % plot the permutations 
                hold on; boundedline(bincent, mean(mean(chanresp_permsR(:,:,34:70),1),3),std(squeeze(chanresp_permsR(:,34:70)))./sqrt(30)*1.96 ,'k','alpha');

                title(strcat('Reconstruction from Alpha Power in Contra-Posterior Electrodes(', group_name{i}, ')'));
                xlabel('Centered Orientation Channel');
                ylabel('Channel Response');
                legend({'error', '-200-0ms','error','30-300ms','error','300-600ms','error','600-800ms'})

                cd(fpath13);
                saveas(gcf, strcat('SUMO_', subject{l}, '_iem_', group_name{i}));
                saveas(gcf, strcat('SUMO_', subject{l}, '_iem_', group_name{i}, '.png'));   
                save(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{l}, '_epoch1_', group_name{i}), 'dstimes', 'stimlabels', 'chanrespsR', 'chanresp_permsR', 'h'); 
            end
        end
    end
end

%% Merging the IEM for group analysis(both left and right AMI)
clearvars -except keep fpath8 & fpath13 & subject & l & group_name;

group_name = {'left AMI', 'right AMI', 'left UMI', 'right UMI'};

cd(fpath13);

for i = 1 : length(subject)
    for j = 1:2
        load(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{i}, '_epoch1_', group_name{j}, '.mat'));
    end
    
    %squeeze the channel response
    tmp1 = mean(chanrespsL,1);
    chanrespsL = squeeze(tmp1);

    tmp2 = mean(chanrespsR,1);
    chanrespsR = squeeze(tmp2);

    tmp3 = mean(chanresp_permsL,1);
    chanresp_permsL = squeeze(tmp3);

    tmp4 = mean(chanresp_permsR,1);
    chanresp_permsR = squeeze(tmp4);

    temp1_CR1L(i,:,:) = chanrespsL;
    temp1_CR1R(i,:,:) = chanrespsR;
    temp1_CRP1L(i,:,:) = chanresp_permsL;
    temp1_CRP1R(i,:,:) = chanresp_permsR;

end

chanrespsL_allsub = temp1_CR1L;
chanrespsR_allsub = temp1_CR1R;
chanresp_permsL_allsub = temp1_CRP1L;
chanresp_permsR_allsub = temp1_CRP1R;

save(strcat('IEM_Exp1_contraposterior_allsubjects_epoch1_ami.mat'), 'dstimes', 'stimlabels', 'chanrespsL_allsub', 'chanresp_permsL_allsub', 'chanrespsR_allsub', 'chanresp_permsR_allsub'); 

%% Plotting Left 
% plotting the group level is different from plotting the single subject
nbins = 7;
binedges = linspace(1,181,nbins+1);
bincent = round(mean([binedges(1:end-1);binedges(2:end)]));

ChanrespsL_allsub = squeeze(nanmean(chanrespsL_allsub,1));

% plot the encoding model on specific time windows with banded lines(left)
% Generally, boundedline(bincent, mean of the timepoint on bins(7 * 1 double), ...
%                           std on the subjects(1st dimension), of the mean on the timepoint(3rd dimension))
tw = dsearchn(dstimes',[-200 -20]'); %baseline, k = gray
figure; boundedline(bincent, squeeze(nanmean(ChanrespsL_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsL_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'c','alpha');

tw = dsearchn(dstimes',[30 300]'); %r = red
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

title('IEM Left Stimulus Reconstruction from Contraposterior Electrodes(TMS1 AMI)(Alpha)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')
legend({'error', '-200-0ms','error','30-300ms','error','300-600ms','error','600-800ms'})

title('IEM Left Stimulus Reconstruction Contraposterior Electrodes(TMS1 AMI)(Alpha)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')

saveas(gcf, 'allsubjects_iem_left_AMI_alpha.png');
%% Plotting Right

ChanrespsR_allsub = squeeze(nanmean(chanrespsR_allsub,1));
%plot the encoding model on specific time windows with banded lines(right)
tw = dsearchn(dstimes',[-200 -20]'); %baseline, k = gray
figure; boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'c','alpha');

tw = dsearchn(dstimes',[30 300]'); %r = red
boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'r','alpha');

tw = dsearchn(dstimes',[300 600]'); %m = magenta
boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'m','alpha');

tw = dsearchn(dstimes',[600 800]'); %b = blue
boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'b','alpha');

% plot the permutations % why is it 34 : 70, and why is it sqrt(30) * 1.96
boundedline(bincent, squeeze(mean(nanmean(chanresp_permsR_allsub(:,:,34:70),3),1)),squeeze(std(mean(chanresp_permsR_allsub(:,:,34:70),3)))./sqrt(30)*1.96,'k','alpha');

hold on;
ylim([-.2 .6]);

title('IEM Right Stimulus Reconstruction from Contraposterior Electrodes(TMS1 AMI)(Alpha)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')
legend({'error', '-200-0ms','error','30-300ms','error','300-600ms','error','600-800ms'})

title('IEM Right Stimulus Reconstruction Contraposterior Electrodes(TMS1 AMI)(Alpha)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')

saveas(gcf, 'allsubjects_iem_right_AMI_alpha.png');

%% Plottingg the heatmaps for the left and right AMI
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

    % plot the dots for the significant timepoints in the heatmap
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
        title('Reconstruction from Alpha Power in ContraPosterior Channels Left AMI(TMS1)', 'Fontsize', titleFontSize)
        xlabel('Time from TMS onset(msec)','Fontsize', axisFontSize)
        ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)
        savefig('allsubjects_iem_left_ami_heatmap');
        saveas(gcf, 'allsubjects_iem_left_ami_heatmap.png');
    elseif k == 2 
        ylabel(h, 'Channel Response (µV^2)','Fontsize', axisFontSize)
        title('Reconstruction from Alpha Power in ContraPosterior Channels Right AMI(TMS1)', 'Fontsize', titleFontSize)
        xlabel('Time from TMS onset(msec)','Fontsize', axisFontSize)
        ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)
        savefig('allsubjects_iem_right_ami_heatmap');
        saveas(gcf, 'allsubjects_iem_right_ami_heatmap.png');
    end
end
%% Merging the IEM for group analysis(both left and right UMI)
clearvars -except keep fpath8 & fpath13 & subject & l & group_name;

group_name = {'left AMI', 'right AMI', 'left UMI', 'right UMI'};

cd(fpath13);

for i = 2 : length(subject)
    for j = 3:4        
        load(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{i}, '_epoch1_', group_name{j}, '.mat'));
    end
    
    %squeeze the channel response
    tmp1 = mean(chanrespsL,1);
    chanrespsL = squeeze(tmp1);

    tmp2 = mean(chanrespsR,1);
    chanrespsR = squeeze(tmp2);

    tmp3 = mean(chanresp_permsL,1);
    chanresp_permsL = squeeze(tmp3);

    tmp4 = mean(chanresp_permsR,1);
    chanresp_permsR = squeeze(tmp4);

    temp1_CR1L(i,:,:) = chanrespsL;
    temp1_CR1R(i,:,:) = chanrespsR;
    temp1_CRP1L(i,:,:) = chanresp_permsL;
    temp1_CRP1R(i,:,:) = chanresp_permsR;

end

chanrespsL_allsub = temp1_CR1L;
chanrespsR_allsub = temp1_CR1R;
chanresp_permsL_allsub = temp1_CRP1L;
chanresp_permsR_allsub = temp1_CRP1R;

save(strcat('IEM_Exp1_contraposterior_allsubjects_epoch1_umi_iem.mat'), 'dstimes', 'stimlabels', 'chanrespsL_allsub', 'chanresp_permsL_allsub', 'chanrespsR_allsub', 'chanresp_permsR_allsub'); 

%% Plotting Left 
% plotting the group level is different from plotting the single subject
nbins = 7;
binedges = linspace(1,181,nbins+1);
bincent = round(mean([binedges(1:end-1);binedges(2:end)]));

ChanrespsL_allsub = squeeze(nanmean(chanrespsL_allsub,1));

% plot the encoding model on specific time windows with banded lines(left)
% Generally, boundedline(bincent, mean of the timepoint on bins(7 * 1 double), ...
%                           std on the subjects(1st dimension), of the mean on the timepoint(3rd dimension))
tw = dsearchn(dstimes',[-200 -20]'); %baseline, k = gray
figure; boundedline(bincent, squeeze(nanmean(ChanrespsL_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsL_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'c','alpha');

tw = dsearchn(dstimes',[30 300]'); %r = red
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

title('IEM Left Stimulus Reconstruction from Contraposterior Electrodes(TMS1 UMI)(Alpha)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response')
legend({'error', '-200-0ms','error','30-300ms','error','300-600ms','error','600-800ms'})

title('IEM Left Stimulus Reconstruction Contraposterior Electrodes(TMS1 UMI)(Alpha)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response')

saveas(gcf, 'allsubjects_iem_left_UMI_alpha.png');
%% Plotting Right

ChanrespsR_allsub = squeeze(nanmean(chanrespsR_allsub,1));
%plot the encoding model on specific time windows with banded lines(right)
tw = dsearchn(dstimes',[-200 -20]'); %baseline, k = gray
figure; boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'c','alpha');

tw = dsearchn(dstimes',[30 300]'); %r = red
boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'r','alpha');

tw = dsearchn(dstimes',[300 600]'); %m = magenta
boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'m','alpha');

tw = dsearchn(dstimes',[600 800]'); %b = blue
boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'b','alpha');

% plot the permutations % why is it 34 : 70, and why is it sqrt(30) * 1.96
boundedline(bincent, squeeze(mean(nanmean(chanresp_permsR_allsub(:,:,34:70),3),1)),squeeze(std(mean(chanresp_permsR_allsub(:,:,34:70),3)))./sqrt(30)*1.96,'k','alpha');

hold on;
ylim([-.2 .6]);

title('IEM Right Stimulus Reconstruction from Contraposterior Electrodes(TMS1 UMI)(Alpha)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response')
legend({'error', '-200-0ms','error','30-300ms','error','300-600ms','error','600-800ms'})

title('IEM Right Stimulus Reconstruction Contraposterior Electrodes(TMS1 UMI)(Alpha)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response')

saveas(gcf, 'allsubjects_iem_right_UMI_alpha.png');

%% Plottingg the heatmaps for the left and right UMI
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

    % plot the dots for the significant timepoints in the heatmap
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
        title('Reconstruction from Alpha Power in ContraPosterior Channels Left UMI(TMS1)', 'Fontsize', titleFontSize)
        xlabel('Time from TMS onset(msec)','Fontsize', axisFontSize)
        ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)
        savefig('allsubjects_iem_left_umi_heatmap');
        saveas(gcf, 'allsubjects_iem_left_umi_heatmap.png');
    elseif k == 2 
        ylabel(h, 'Channel Response (µV^2)','Fontsize', axisFontSize)
        title('Reconstruction from Alpha Power in ContraPosterior Channels Right UMI(TMS1)', 'Fontsize', titleFontSize)
        xlabel('Time from TMS onset(msec)','Fontsize', axisFontSize)
        ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)
        savefig('allsubjects_iem_right_umi_heatmap');
        saveas(gcf, 'allsubjects_iem_right_umi_heatmap.png');
    end
end

% finished
load handel;
sound(y, Fs);