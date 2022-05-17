clear;
close all;

%% Start the Journey
%%%%%%%%%%%%THIS ONE DOES NOT USE THE SUPERTRIALS BUT ORIGINAL EEG DATA%%%%%%%%%%%%%%%%%%%%%%%

% Generally, there are five parts for the EEG processing
% 1. eeg_preproc_part1_merging -> merging beh data onto the eeg raw data
% 2. eeg_preproc_part2_Epoching -> loading channel locations and epoching eeg for analysis
%
% Start to run locally and manually
% 3. eeg_preproc_part3_Demeaning -> calculating pre-stim baseline and removing DC offset
% 4. eeg_preproc_part4_prepareICA -> preprocess eeg data for later fast_ICA
% 5. eeg_preproc_part5_ICA -> running on fast_ICA twice

% 6. eeg_analysis_n2pc -> calculating n2pc on 28 channels
% 7. eeg_analysis_iem_voltage -> iem on probe/tms/cue/stim data (from -300ms to 900ms)
% 8. eeg_analysis_iem_a_power -> iem on probe/tms/cue/qstat stim data (from -300ms to 900ms)

% In this part, I will run the iem on Alpha-power in probe1 through processed
% data

% This script requires EEGLAB v2020.0 in MATLAB

% The flow for this part:
%                                Loading EEG sets
%                                       |
%                              Prepare for the IEM
%                                       |
%                                 IEM on Alpha-power
%                                       |
%                                 Group Analysis  
%                                       |
%                                    Save out

%% Load the Basics

addpath(genpath('E:\SUMO_further_data_pack_zx\ERP_TESA\envi'));
fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';
fpath13 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\test';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120', ...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};
        
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

group_name = {'left stim', 'right stim'};

%% Start the IEM
for l = 1:length(subject)
    %1 : length(subject)
    %for t = 1 : length(type) % Picking different types
    for t = 1 % This time we only preprocess the Probe1 data
        %% Step 1: Load EEG data
        clearvars -except keep fpath8 & fpath13 & subject & l & type & group_name & t;
        
        cd(fpath8);
        EEG = pop_loadset('filename', strcat(subject{l}, '_before_iem_', type{t}, '.set'), 'filepath', fpath8);
                
        EEG_data_doubled = double(EEG.data);% change from single dataset to double
        p = struct2table(EEG.event); % change the EEG.event to a cell array and then categorize into different bins 
        
        
        %% Step 3: filter trials into two groups: left_cue, right_cue
        %get the event trial number; target location1; test orientation; epoch#; left orientation(stimulus); right orientaion(stimulus); 
        eventBins = [p.trialnum, p.targetlocation, p.targetorient, p.epoch, p.leftori, p.rightori];
        
        %Filter trials into left and right group
        %bins(1.trial#, 2.trialloc, 3.targetori, 4.epoch#, 5.leftori, 6.rightori)
        
        nbins = 7;
        binedges = round(linspace(1,181,nbins+1));
        bincent = round(mean([binedges(1:end-1);binedges(2:end)]));
        [leftTrial, leftStimlabels, rightTrial, rightStimlabels]  = filter_bins_for_iem_stim(p, binedges, bincent);
        
        %% Step 4: IEM
        groups = {leftTrial, rightTrial};
        labels = {leftStimlabels, rightStimlabels};
        alpha_power = zeros(size(EEG.data)); % preallocate the alpha_power dataset
        pow_base = [];

        % Run the dothewave function to have pow matrix
        % We can change from [2 50] to [8 13](alpha power)
        % This time we pick the pre-TMS interval as baseline for FFT 
        [pow, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [8 13], 6, 4, 1, [-200 -20], EEG.times);

        % Average the 6 frequencies into one
        alpha_power = squeeze(mean(pow,2));   
        
        for i = 1 : length(group_name) % notice that we used the parallel envi in permutation
            chanrespsL = [];
            chanrespsR = [];
            chanresp_permsL = [];
            chanresp_permsR = [];
            
            clear trials & stimlabels;
            
            [super_charge, stimlabels, h] = supertrial_3trial(groups{i}, alpha_power);
            
            clear chanresp & weights & dstimes & chanresp_perm;
            
            % for all posterior electrodes
            % impchan = [23,53,22,54,21,51,19,50,20,48,49,18,10,41,11,42,12,15,44,14,43,45,46,16]; 

            %for contraposterior electrodes
            switch i
                case 1
                    impchan = [23,53,22,54,21,51,19,50,20,48,49,18]; %channels in R hem
                otherwise
                    impchan = [10,41,11,42,12,15,44,14,43,45,46,16]; %channels in L hem
            end
        
            [chanresp, ~, dstimes] = iemori(super_charge(impchan,:,:),stimlabels,4,EEG.times);

            % permutation test to get null hypothesis - takes a LONG time
            % in the permutation, the stimlabels need to be shuffled
            nperm = 100;
            parfor p = 1:nperm
                disp(p)
                [tmp, ~, dstimes] = iemori(super_charge(impchan,:,:),stimlabels(randperm(length(stimlabels))),4,EEG.times); 
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
                if h(1,1,3) == 0  
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
                tw = dsearchn(dstimes',[-80 0]'); %baseline, k = gray
                figure; boundedline(bincent, nanmean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2)))./sqrt(7)),'c','alpha'); 

                tw = dsearchn(dstimes',[20 200]'); %r = red
                hold on; boundedline(bincent, mean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2))))./sqrt(7),'r','alpha');

                tw = dsearchn(dstimes',[200 400]'); %m = magenta
                hold on; boundedline(bincent, mean(ChanrespsL(:,tw(1):tw(2)),2),std(squeeze(chanrespsL(:,tw(1):tw(2))))./sqrt(7) ,'m','alpha');
                
                % plot the permutations 
                hold on; boundedline(bincent, mean(mean(chanresp_permsL(:,:,34:70),1),3),std(squeeze(chanresp_permsL(:,34:70)))./sqrt(30)*1.96 ,'k','alpha');
                
                hold on;
                ylim([-.2 .6]);
                
                title(strcat('Reconstruction from Alpha-power in Contra-Posterior Electrodes(', group_name{i}, ')'));
                xlabel('Centered Orientation Channel');
                ylabel('Channel Response');
                legend({'error', '-80-0ms','error','20-200ms','error','200-400ms'})

                cd(fpath13);
                saveas(gcf, strcat('SUMO_', subject{l}, '_iem_', group_name{i}, '.png'));    
                save(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{l}, '_epoch1_', group_name{i}), 'dstimes', 'stimlabels', 'chanrespsL', 'chanresp_permsL'); 
            else % for right stimuli
                avgacrosstrials = mean(chanresp,3); %create a variable that is now 2D
                avgacrosstrials_perm = mean(chanresp_perm,3); %create a variable that is now 2D
            
                % checking if the subject miss the #3 channel in EEG.data
                % for subjects who do not have the 3rd channel, we copy the 5th channel to the 3rd
                if h(1,1,3) == 0  
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
                tw = dsearchn(dstimes',[-80 0]'); %baseline, k = gray
                figure; boundedline(bincent, nanmean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2)))./sqrt(7)),'c','alpha'); 

                tw = dsearchn(dstimes',[20 200]'); %r = red
                hold on; boundedline(bincent, mean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2))))./sqrt(7),'r','alpha');

                tw = dsearchn(dstimes',[200 400]'); %m = magenta
                hold on; boundedline(bincent, mean(ChanrespsR(:,tw(1):tw(2)),2),std(squeeze(chanrespsR(:,tw(1):tw(2))))./sqrt(7) ,'m','alpha');

                % plot the permutations 
                hold on; boundedline(bincent, mean(mean(chanresp_permsR(:,:,34:70),1),3),std(squeeze(chanresp_permsR(:,34:70)))./sqrt(30)*1.96 ,'k','alpha');
                
                hold on;
                ylim([-.2 .6]);
                
                title(strcat('Reconstruction from Alpha-power in Contra-Posterior Electrodes(', group_name{i}, ')'));
                xlabel('Centered Orientation Channel');
                ylabel('Channel Response');
                legend({'error', '-80-0ms','error','20-200ms','error','200-400ms'})

                cd(fpath13);
                saveas(gcf, strcat('SUMO_', subject{l}, '_iem_', group_name{i}, '.png'));    
                save(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{l}, '_epoch1_', group_name{i}), 'dstimes', 'stimlabels', 'chanrespsR', 'chanresp_permsR'); 
            end
        end
    end
end

%% Merging the IEM for group analysis(both left and right cue)
clearvars -except keep fpath8 & fpath13 & subject & l & group_name;

group_name = {'left stim', 'right stim'};

cd(fpath13);

for i = 1 : length(subject)
    for j = 1:2
        load(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{i}, '_epoch1_', group_name{j}, '.mat'));
    end
    
    %squeeze the channel response
    chanrespsL = squeeze(mean(chanrespsL,1));
    chanrespsR = squeeze(mean(chanrespsR,1));
    chanresp_permsL = squeeze(mean(chanresp_permsL,1));
    chanresp_permsR = squeeze(mean(chanresp_permsR,1));

    temp1_CR1L(i,:,:) = chanrespsL;
    temp1_CR1R(i,:,:) = chanrespsR;
    temp1_CRP1L(i,:,:) = chanresp_permsL;
    temp1_CRP1R(i,:,:) = chanresp_permsR;

end

chanrespsL_allsub = temp1_CR1L;
chanrespsR_allsub = temp1_CR1R;
chanresp_permsL_allsub = temp1_CRP1L;
chanresp_permsR_allsub = temp1_CRP1R;

save(strcat('IEM_Exp1_contraposterior_allsubjects_epoch1_stim_iem.mat'), 'dstimes', 'stimlabels', 'chanrespsL_allsub', 'chanresp_permsL_allsub', 'chanrespsR_allsub', 'chanresp_permsR_allsub'); 

%% Plotting Left 
% plotting the group level is different from plotting the single subject
nbins = 7;
binedges = linspace(1,181,nbins+1);
bincent = round(mean([binedges(1:end-1);binedges(2:end)]));

ChanrespsL_allsub = squeeze(nanmean(chanrespsL_allsub,1));

% plot the encoding model on specific time windows with banded lines(left)
% Generally, boundedline(bincent, mean of the timepoint on bins(7 * 1 double), ...
%                           std on the subjects(1st dimension), of the mean on the timepoint(3rd dimension))
tw = dsearchn(dstimes',[-80 0]'); %baseline, k = gray
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
% ylim([-.2 .6]);

title('IEM Left Stimulus Reconstruction from Contraposterior Electrodes(Stim)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')
legend({'error', '-80-0ms','error','20-200ms','error','200-400ms'})

title('IEM Left Stimulus Reconstruction Contraposterior Electrodes(Stim)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')

saveas(gcf, 'allsubjects_iem_left_stim.png');
%% Plotting Right

ChanrespsR_allsub = squeeze(nanmean(chanrespsR_allsub,1));
%plot the encoding model on specific time windows with banded lines(right)
tw = dsearchn(dstimes',[-80 0]'); %baseline, k = gray
figure; boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'c','alpha');

tw = dsearchn(dstimes',[20 200]'); %r = red
boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'r','alpha');

tw = dsearchn(dstimes',[200 400]'); %m = magenta
boundedline(bincent, squeeze(nanmean(ChanrespsR_allsub(:,tw(1):tw(2)),2)),squeeze(std(mean(chanrespsR_allsub(:,:,tw(1):tw(2)),3))./sqrt(10)),'m','alpha');

% plot the permutations % why is it 34 : 70, and why is it sqrt(30) * 1.96
boundedline(bincent, squeeze(mean(nanmean(chanresp_permsR_allsub(:,:,34:70),3),1)),squeeze(std(mean(chanresp_permsR_allsub(:,:,34:70),3)))./sqrt(30)*1.96,'k','alpha');

hold on;
% ylim([-.2 .6]);

title('IEM Right Stimulus Reconstruction from Contraposterior Electrodes(Stim)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')
legend({'error', '-80-0ms','error','20-200ms','error','200-400ms'})

title('IEM Right Stimulus Reconstruction Contraposterior Electrodes(Stim)')
xlabel('Centered Orientation Channel (°)')
ylabel('Channel Response (µV^2)')

saveas(gcf, 'allsubjects_iem_right_stim.png');

%% Plottingg the heatmaps for the left and right cue
% Load Parameters
centerind = 4;
nchan = 7;

% Set up title font size
titleFontSize = 20;
axisFontSize = 14;
textFontSize = 12;

% Smooth the data to 50ms window by Gaussian filter
for k =  1 : 2
    if k == 1
%         n = size(chanrespsL_allsub);
%         m = n(1); % # of bins(6 bins or 7)
%         n = n(2); % # of time points
%         % Smoothing data per channel per trial, accross the 1200ms interval
%         for i = 1 : m % # of subjects
%             for j = 1 : n % # of bins
%                 % The channel response is 1200ms (300 time points)
%                 % Use 6 as the window to get 50ms per sample
%                 chanrespsL_allsub(i,j,:) = smoothdata(chanrespsL_allsub(i,j,:), 'gaussian', 6);
%             end
%         end
        cond = chanrespsL_allsub;
    else
%         n = size(chanrespsR_allsub);
%         m = n(1); % # of bins(6 bins or 7)
%         n = n(2); % # of time points
%         % Smoothing data per channel per trial, accross the 1200ms interval
%         for i = 1 : m % # of subjects
%             for j = 1 : n % # of bins
%                 chanrespsR_allsub(i,j,:) = smoothdata(chanrespsR_allsub(i,j,:), 'gaussian', 6);
%             end
%         end
        cond = chanrespsR_allsub;
    end

    % the process for generating the heat map
    orioffset = linspace(-90, 90, nchan); %for plotting
    figure; contourf(dstimes, orioffset, squeeze(mean(cond, 1)),30,'linec','none'); %make a contour plot
    xlim([min(dstimes) max(dstimes)]); %can trim the time range of the plot here if you want. Right now setting to min and max does nothing.
    h = colorbar;
    caxis([-0.1 0.6]); %scale the colorbar axis to just region of interest (-0.1 to 0.6)

    % Opening the figure in the full size
    set(gcf, 'Position', get(0, 'Screensize'));
    % determine if it is left or right channel response
    if k == 1
        ylabel(h, 'Channel Response (µV^2)')
        title('Reconstruction from Alpha-power in ContraPosterior Channels Left Stimulus(Stim)', 'Fontsize', titleFontSize)
        xlabel('Time from cue onset(msec)','Fontsize', axisFontSize)
        ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)
        savefig('allsubjects_iem_left_stim_heatmap');
        saveas(gcf, 'allsubjects_iem_left_stim_heatmap.png');
    elseif k == 2 
        ylabel(h, 'Channel Response (µV^2)')
        title('Reconstruction from Alpha-power in ContraPosterior Channels Right Stimulus(Stim)', 'Fontsize', titleFontSize)
        xlabel('Time from cue onset(msec)','Fontsize', axisFontSize)
        ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)
        savefig('allsubjects_iem_right_stim_heatmap');
        saveas(gcf, 'allsubjects_iem_right_stim_heatmap.png');
    end
end
% finished