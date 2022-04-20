function [super_charge, stimlabels, h] = supertrial_3trial(bin, EEG_data)

%Input : bin : trial numbers of AMI and UMI trials
%        EEG_data: raw eeg_data, has been doubled before

%Output: super_charge: eeg data have been categorized into 
%                      different bins based on target orientation
%        stimlabels: matrices of stimlabels (14 degrees to 181 degrees)
%                    for later iem
%        


super_trial_index = 1;

%get the trial number for each bins

% a is the number of trials
% h is a 3-d matrix 
% h[index number(from 1 to ...), trial number from EEG(p.epoch), bin number(1:7)]
for i = 1 : 7
    a = 0;
    for j = 1 : length(bin)
        if bin(j,i) ~= 0
            a = a + 1;
            h(a,1,i) = a;
            h(a,2,i) = bin(j,i);
        end
    end
    
    % b is the number of supertrials
    b = fix(a/3); 
    c = mod(a,3);

    % get the super trials into a 4-D matrix (channels x time x supertrial number x bin)
    % for instance, if there are 13 trials for bin 1, there will be 4 supertrials: (ave(trial 1,5,9), ave(2,6,10), ave(3,7,11), ave(4,8,12,13))
    % similarly, if there are 14 trilas for bin 1, there will be 4 supertrials: (ave(trial 1,5,9), ave(2,6,10), ave(3,7,11), ave(4,8,12,13,14))
    for k = 1 : b
        if k < b
            B(:,:,k,i) = mean(EEG_data(:,:,[h(k,2,i),h(k+b,2,i),h(k+b*2,2,i)]),3);
        elseif k == b
            switch c
                case 0
                    B(:,:,k,i) = mean(EEG_data(:,:,[h(k,2,i),h(k+b,2,i),h(k+b*2,2,i)]),3);
                case 1
                    B(:,:,k,i) = mean(EEG_data(:,:,[h(k,2,i),h(k+b,2,i),h(k+b*2,2,i),h(k+b*2+1,2,i)]),3);
                otherwise
                %here are some changes 
                    B(:,:,k,i) = mean(EEG_data(:,:,[h(k,2,i),h(k+b,2,i),h(k+b*2,2,i),h(k+b*2+1,2,i),h(k+b*2+2,2,i)]),3);
            end
        end
        stimlabels(1,super_trial_index) = i;%later will be changed to the center degree of bins
        super_trial_index = super_trial_index + 1; %the number of supertrials in one bin
    end
end

% get the change the bin number to bin center orientation and then get the stimlabels for later IEM
for i = 1 : length(stimlabels)
    if stimlabels(1,i) == 1
        stimlabels(1,i) = 14;
    elseif stimlabels(1,i) == 2
        stimlabels(1,i) = 40;
    elseif stimlabels(1,i) == 3
        stimlabels(1,i) = 65;
    elseif stimlabels(1,i) == 4
        stimlabels(1,i) = 91;
    elseif stimlabels(1,i) == 5
        stimlabels(1,i) = 117;
    elseif stimlabels(1,i) == 6
        stimlabels(1,i) = 142;
    else
        stimlabels(1,i) = 168;
    end
end

d = size(B);
g = 1;

% extract the supertrials eeg data from original EEG.data
for i = 1 : 7
    for m = 1 : d(3)
        if B(1,1,m,i) ~= 0
            super_charge(:,:,g) = B(:,:,m,i);
            g = g + 1;
        end
    end
end