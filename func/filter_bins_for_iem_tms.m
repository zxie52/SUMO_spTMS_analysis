function [leftAMI, leftUMI, rightAMI, rightUMI]  = filter_bins_for_iem_tms(p, binedges, bincent)

% The input is the struct2table(EEG.event), binedges, bincent
% The output is trial number and bin# for left/right AMI and UMI


for i =  1 : height(p)
    j = 1;
    drop_case = [];
    % filter the left stimulus orientations into seven bins 
    switch p.leftori(i)
        % use round functions to have the integer
        case num2cell(round(binedges(1):binedges(2)))
            p.stimlabels(i) = bincent(1);
            p.leftBin(i) = 1;
        case num2cell(round(binedges(2):binedges(3)))
            p.stimlabels(i) = bincent(2);
            p.leftBin(i) = 2;
        case num2cell(round(binedges(3):binedges(4)))
            p.stimlabels(i) = bincent(3);
            p.leftBin(i) = 3;
        case num2cell(round(binedges(4):binedges(5)))
            p.stimlabels(i) = bincent(4);
            p.leftBin(i) = 4;
        case num2cell(round(binedges(5):binedges(6)))
            p.stimlabels(i) = bincent(5);
            p.leftBin(i) = 5;
        case num2cell(round(binedges(6):binedges(7)))
            p.stimlabels(i) = bincent(6);
            p.leftBin(i) = 6;
        case num2cell(round(binedges(7):binedges(8)))
            p.stimlabels(i) = bincent(7);
            p.leftBin(i) = 7;
        otherwise 
            fprintf("Something is wrong on the orientations\n");
            fprintf("The wrong trial is %d\n", p.epoch(i));
            drop_case(j) = p.epoch(i);
            j = j + 1;     
    end
    % Then filter the right stimulus orientation to seven bins
    switch p.rightori(i)
        case num2cell(round(binedges(1):binedges(2)))
            p.stimlabels(i) = bincent(1);
            p.rightBin(i) = 1;
        case num2cell(round(binedges(2):binedges(3)))
            p.stimlabels(i) = bincent(2);
            p.rightBin(i) = 2;
        case num2cell(round(binedges(3):binedges(4)))
            p.stimlabels(i) = bincent(3);
            p.rightBin(i) = 3;
        case num2cell(round(binedges(4):binedges(5)))
            p.stimlabels(i) = bincent(4);
            p.rightBin(i) = 4;
        case num2cell(round(binedges(5):binedges(6)))
            p.stimlabels(i) = bincent(5);
            p.rightBin(i) = 5;
        case num2cell(round(binedges(6):binedges(7)))
            p.stimlabels(i) = bincent(6);
            p.rightBin(i) = 6;
        case num2cell(round(binedges(7):binedges(8)))
            p.stimlabels(i) = bincent(7);
            p.rightBin(i) = 7;
        otherwise 
            fprintf("Something is wrong on the orientations\n");
            fprintf("The wrong trial is %d\n", p.epoch(i));
            drop_case(j) = p.epoch(i);
            j = j + 1;     
    end
end

p(drop_case, :) = [];
fprintf("The wrong trial has been dropped! Filtering has been finished\n");

% separate different conditions and preallocate the space
% ntrials on left AMI should be same with those on right UMI and vice versa
leftAMI = zeros(length(find(p.targetlocation == 1)), 7);
leftUMI = zeros(length(find(p.targetlocation == 2)), 7);
rightAMI = zeros(length(find(p.targetlocation == 2)), 7);
rightUMI = zeros(length(find(p.targetlocation == 1)), 7);

% % identify if this trial is cue1 or cue2
% if strcmp(type, 'cue1')
%     loc = p.targetlocation;
% elseif strcmp(type, 'cue2')
%     loc = p.targetlocation2;
% end
loc = p.targetlocation;
% create the 2-d bins matrix: ntrials x nbins
for i = 1 : height(p)
    switch loc(i)
        case 1
            leftAMI(i, p.leftBin(i)) = p.epoch(i);
            rightUMI(i, p.rightBin(i)) = p.epoch(i);
        case 2
            rightAMI(i, p.rightBin(i)) = p.epoch(i);
            leftUMI(i, p.leftBin(i)) = p.epoch(i);
        otherwise
            fprintf("Something went wrong");
    end
end

% The output would be 2-d matrix, with the cell bin the trial number for specific bin
% for example, if the leftAMI(2,5) = 2, then the 2nd trial is leftAMI and bin should be #5 bin

