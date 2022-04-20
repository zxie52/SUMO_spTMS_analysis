function [leftTrial, leftStimlabels, rightTrial, rightStimlabels]  = filter_bins_for_iem_probe(str, cond, p, binedges, bincent)

% The input is the struct2table(EEG.event), binedges, bincent
% The output is the left and right trial numbers and compatible stimlabels

% identify if this trial is cue1 or cue2
if strcmp(cond, 'probe1')
    loc = p.targetlocation;
    if strcmp(str, 'testorient')
        orient = p.testorient;
    elseif strcmp(str, 'targetorient')
        orient = p.targetorient;
    end
elseif strcmp(cond, 'probe2')
    loc = p.targetlocation2;
    if strcmp(str, 'testorient')
        orient = p.testorient2;
    elseif strcmp(str, 'targetorient')
        orient = p.targetorient2;
    end
end

for i =  1 : height(p)
    j = 1;
    drop_case = [];
    switch orient(i)
        case num2cell(round(binedges(1)-1:binedges(2)))
            p.stimlabels(i) = bincent(1);
            p.bin(i) = 1;
        case num2cell(round(binedges(2):binedges(3)))
            p.stimlabels(i) = bincent(2);
            p.bin(i) = 2;
        case num2cell(round(binedges(3):binedges(4)))
            p.stimlabels(i) = bincent(3);
            p.bin(i) = 3;
        case num2cell(round(binedges(4):binedges(5)))
            p.stimlabels(i) = bincent(4);
            p.bin(i) = 4;
        case num2cell(round(binedges(5):binedges(6)))
            p.stimlabels(i) = bincent(5);
            p.bin(i) = 5;
        case num2cell(round(binedges(6):binedges(7)))
            p.stimlabels(i) = bincent(6);
            p.bin(i) = 6;
        case num2cell(round(binedges(7):binedges(8)))
            p.stimlabels(i) = bincent(7);
            p.bin(i) = 7;
        otherwise 
            fprintf("Something is wrong on the orientations\n");
            fprintf("The wrong trial is %d\n", p.epoch(i));
            drop_case(j) = p.epoch(i);
            j = j + 1;
    end
end

p(drop_case, :) = [];
fprintf("The wrong trial is dropped!\n");

p = movevars(p, 'stimlabels', 'After', 'epoch');

% separate different conditionsn and preallocate the space
leftTrial = zeros(length(p.targetlocation), 7);
leftStimlabels = zeros(length(find(p.targetlocation == 1)), 1);
rightTrial = zeros(length(p.targetlocation), 7);
rightStimlabels = zeros(length(find(p.targetlocation == 2)), 1);

k = 1;
j = 1;
for i = 1 : height(p)
    switch p.targetlocation(i)
        case 1
            leftTrial(i,p.bin(i)) = p.epoch(i);
            leftStimlabels(j) = p.stimlabels(i);
            j = j + 1;
        case 2
            rightTrial(i,p.bin(i)) = p.epoch(i);
            rightStimlabels(k) = p.stimlabels(i);
            k = k + 1;
        otherwise
            fprintf("Something went wrong");
    end
end

