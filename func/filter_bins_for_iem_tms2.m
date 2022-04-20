%%% The purpose of this function is to separate different trials into three
%%% dimensional matrices: 
%%% left_bin_AMI_stay; left_bin_AMI_switch;
%%% left_bin_UMI_stay; left_bin_UMI_switch; right_bin_AMI_stay;
%%% right_bin_AMI_switch; right_bin_UMI_stay; right_bin_UMI_switch

%%% Generally, we first check if the trial is switch or stay;
%%% Then, check if it is left or right
%%% Same with the TMS1, leftAMI and rightUMI will be in the same trial number, 
%%% vice versa

%%% The input bins should have seven columns:
%%% column1: trial number;
%%% column2: target location 1;
%%% column3: target location 2;
%%% column4: target orientation 2;
%%% column5: epoch #;
%%% column6: left orientation (stimulus S 48)
%%% column7: right orientation (stimulus S 48)

%%% Input: EEG.event, bindeges and bincent

%%% Output: eight matrices, which will be used in later supertrial function

%%% This function is still a little bit redundant, should be refined

function [left_bin_AMI_stay, left_bin_AMI_switch, ...
          left_bin_UMI_stay, left_bin_UMI_switch, ...
          right_bin_AMI_stay, right_bin_AMI_switch, ...
          right_bin_UMI_stay, right_bin_UMI_switch] = filter_bins_for_iem_tms2(p, binedges)

% preallocate the space for memory
left_bin_AMI_stay = zeros(height(p),7);
left_bin_AMI_switch = zeros(height(p),7);
left_bin_UMI_stay = zeros(height(p),7);
left_bin_UMI_switch = zeros(height(p),7);

right_bin_AMI_stay = zeros(height(p),7);
right_bin_AMI_switch = zeros(height(p),7);
right_bin_UMI_stay = zeros(height(p),7);
right_bin_UMI_switch = zeros(height(p),7);

for i = 1 : height(p)
    if p.targetlocation(i) == p.targetlocation2(i) % this trial should be stay
        % stay trials
        if p.targetlocation2(i) == 1 
            %% this trial is left: leftAMI stay and rightUMI stay
            if ismember(p.targetorient2(i), round(binedges(1):binedges(2)))
                left_bin_AMI_stay(i,1) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(2):binedges(3)))
                left_bin_AMI_stay(i,2) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(3):binedges(4)))
                left_bin_AMI_stay(i,3) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(4):binedges(5)))
                left_bin_AMI_stay(i,4) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(5):binedges(6)))
                left_bin_AMI_stay(i,5) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(6):binedges(7)))
                left_bin_AMI_stay(i,6) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(7):binedges(8)))
                left_bin_AMI_stay(i,7) = p.epoch(i);
            end
            
            %% for the right UMI stay(right ori)
            if ismember(p.rightori(i), round(binedges(1):binedges(2)))
               right_bin_UMI_stay(i,1) = p.epoch(i); 
            end
            if ismember(p.rightori(i), round(binedges(2):binedges(3)))
                right_bin_UMI_stay(i,2) = p.epoch(i);
            end
            if ismember(p.rightori(i), round(binedges(3):binedges(4)))
                right_bin_UMI_stay(i,3) = p.epoch(i);
            end
            if ismember(p.rightori(i), round(binedges(4):binedges(5)))
                right_bin_UMI_stay(i,4) = p.epoch(i);
            end
            if ismember(p.rightori(i), round(binedges(5):binedges(6)))
                right_bin_UMI_stay(i,5) = p.epoch(i);
            end
            if ismember(p.rightori(i), round(binedges(6):binedges(7)))
                right_bin_UMI_stay(i,6) = p.epoch(i);
            end
            if ismember(p.rightori(i), round(binedges(7):binedges(8)))
                right_bin_UMI_stay(i,7) = p.epoch(i);
            end 
                
        elseif p.targetlocation2(i) == 2
        %% then there should be rightAMI and leftUMI (for stay groups)
            if ismember(p.targetorient2(i), round(binedges(1):binedges(2)))
                right_bin_AMI_stay(i,1) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(2):binedges(3)))
                right_bin_AMI_stay(i,2) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(3):binedges(4)))
                right_bin_AMI_stay(i,3) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(4):binedges(5)))
                right_bin_AMI_stay(i,4) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(5):binedges(6)))
                right_bin_AMI_stay(i,5) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(6):binedges(7)))
                right_bin_AMI_stay(i,6) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(7):binedges(8)))
                right_bin_AMI_stay(i,7) = p.epoch(i);
            end
        %% then for left_bin_umi(leftori)
            if ismember(p.leftori(i), round(binedges(1):binedges(2)))
                left_bin_UMI_stay(i,1) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(2):binedges(3)))
                left_bin_UMI_stay(i,2) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(3):binedges(4)))
                left_bin_UMI_stay(i,3) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(4):binedges(5)))
                left_bin_UMI_stay(i,4) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(5):binedges(6)))
                left_bin_UMI_stay(i,5) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(6):binedges(7)))
                left_bin_UMI_stay(i,6) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(7):binedges(8)))
                left_bin_UMI_stay(i,7) = p.epoch(i);
            end 
        end
    %% this trial should then be switch    
    elseif p.targetlocation(i) ~= p.targetlocation2(i) 
        if p.targetlocation2(i) == 1 
            %% this trial is left: leftAMI switch and rightUMI switch
            if ismember(p.targetorient2(i), round(binedges(1):binedges(2)))
                left_bin_AMI_switch(i,1) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(2):binedges(3)))
                left_bin_AMI_switch(i,2) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(3):binedges(4)))
                left_bin_AMI_switch(i,3) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(4):binedges(5)))
                left_bin_AMI_switch(i,4) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(5):binedges(6)))
                left_bin_AMI_switch(i,5) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(6):binedges(7)))
                left_bin_AMI_switch(i,6) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(7):binedges(8)))
                left_bin_AMI_switch(i,7) = p.epoch(i);
            end
            %% for the right UMI switch(rightori)
            if ismember(p.rightori(i), round(binedges(1):binedges(2)))
               right_bin_UMI_switch(i,1) = p.epoch(i); 
            end
            if ismember(p.rightori(i), round(binedges(2):binedges(3)))
                right_bin_UMI_switch(i,2) = p.epoch(i);
            end
            if ismember(p.rightori(i), round(binedges(3):binedges(4)))
                right_bin_UMI_switch(i,3) = p.epoch(i);
            end
            if ismember(p.rightori(i), round(binedges(4):binedges(5)))
                right_bin_UMI_switch(i,4) = p.epoch(i);
            end
            if ismember(p.rightori(i), round(binedges(5):binedges(6)))
                right_bin_UMI_switch(i,5) = p.epoch(i);
            end
            if ismember(p.rightori(i), round(binedges(6):binedges(7)))
                right_bin_UMI_switch(i,6) = p.epoch(i);
            end
            if ismember(p.rightori(i), round(binedges(7):binedges(8)))
                right_bin_UMI_switch(i,7) = p.epoch(i);
            end 
            
        elseif p.targetlocation2(i) == 2
        %% then there should be rightAMI and leftUMI (for switch groups)
            if ismember(p.targetorient2(i), round(binedges(1):binedges(2)))
                right_bin_AMI_switch(i,1) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(2):binedges(3)))
                right_bin_AMI_switch(i,2) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(3):binedges(4)))
                right_bin_AMI_switch(i,3) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(4):binedges(5)))
                right_bin_AMI_switch(i,4) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(5):binedges(6)))
                right_bin_AMI_switch(i,5) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(6):binedges(7)))
                right_bin_AMI_switch(i,6) = p.epoch(i);
            end
            if ismember(p.targetorient2(i), round(binedges(7):binedges(8)))
                right_bin_AMI_switch(i,7) = p.epoch(i);
            end
        %% then for left_bin_umi(leftori)
            if ismember(p.leftori(i), round(binedges(1):binedges(2)))
                left_bin_UMI_switch(i,1) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(2):binedges(3)))
                left_bin_UMI_switch(i,2) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(3):binedges(4)))
                left_bin_UMI_switch(i,3) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(4):binedges(5)))
                left_bin_UMI_switch(i,4) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(5):binedges(6)))
                left_bin_UMI_switch(i,5) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(6):binedges(7)))
                left_bin_UMI_switch(i,6) = p.epoch(i);
            end
            if ismember(p.leftori(i), round(binedges(7):binedges(8)))
                left_bin_UMI_switch(i,7) = p.epoch(i);
            end 
        end
    end
end
