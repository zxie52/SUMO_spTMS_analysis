function [leftTrial, leftStimlabels, rightTrial, rightStimlabels]  = filter_bins_for_iem_stim(p, binedges, bincent)

% The input is the struct2table(EEG.event), binedges, bincent
% The output is the left and right trial numbers and compatible stimlabels


for i =  1 : height(p)
    j = 1;
    drop_case = [];
    switch p.leftori(i)
        case num2cell(round(binedges(1):binedges(2)))
            p.stimlabels(i) = bincent(1);
            p.leftoribin(i) = 1;
        case num2cell(round(binedges(2):binedges(3)))
            p.stimlabels(i) = bincent(2);
            p.leftoribin(i) = 2;
        case num2cell(round(binedges(3):binedges(4)))
            p.stimlabels(i) = bincent(3);
            p.leftoribin(i) = 3;
        case num2cell(round(binedges(4):binedges(5)))
            p.stimlabels(i) = bincent(4);
            p.leftoribin(i) = 4;
        case num2cell(round(binedges(5):binedges(6)))
            p.stimlabels(i) = bincent(5);
            p.leftoribin(i) = 5;
        case num2cell(round(binedges(6):binedges(7)))
            p.stimlabels(i) = bincent(6);
            p.leftoribin(i) = 6;
        case num2cell(round(binedges(7):binedges(8)))
            p.stimlabels(i) = bincent(7);
            p.leftoribin(i) = 7;
        otherwise 
            fprintf("Something is wrong on the orientations\n");
            fprintf("The wrong trial is %d\n", p.epoch(i));
            drop_case(j) = p.epoch(i);
            j = j + 1;
    end
    switch p.rightori(i)
        case num2cell(round(binedges(1):binedges(2)))
            p.stimlabels(i) = bincent(1);
            p.rightoribin(i) = 1;
        case num2cell(round(binedges(2):binedges(3)))
            p.stimlabels(i) = bincent(2);
            p.rightoribin(i) = 2;
        case num2cell(round(binedges(3):binedges(4)))
            p.stimlabels(i) = bincent(3);
            p.rightoribin(i) = 3;
        case num2cell(round(binedges(4):binedges(5)))
            p.stimlabels(i) = bincent(4);
            p.rightoribin(i) = 4;
        case num2cell(round(binedges(5):binedges(6)))
            p.stimlabels(i) = bincent(5);
            p.rightoribin(i) = 5;
        case num2cell(round(binedges(6):binedges(7)))
            p.stimlabels(i) = bincent(6);
            p.rightoribin(i) = 6;
        case num2cell(round(binedges(7):binedges(8)))
            p.stimlabels(i) = bincent(7);
            p.rightoribin(i) = 7;
        otherwise 
            fprintf("Something is wrong on the orientations\n");
            fprintf("The wrong trial is %d\n", p.epoch(i));
            drop_case(j) = p.epoch(i);
            j = j + 1;
    end
end

p(drop_case, :) = [];
fprintf("The wrong trial has been dropped! Filtering has been finished\n");

p = movevars(p, 'stimlabels', 'After', 'epoch');

% separate different conditionsn and preallocate the space
leftTrial = zeros(height(p), 7);
leftStimlabels = zeros(height(p), 1);
rightTrial = zeros(height(p), 7);
rightStimlabels = zeros(height(p), 1);

k = 1;
j = 1;
for i = 1 : height(p)
    leftTrial(i,p.leftoribin(i)) = p.epoch(i);
    rightTrial(i,p.rightoribin(i)) = p.epoch(i);
end
