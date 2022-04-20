% Input: IEM channel response and relative permutation

% Output: slope1: the slope changing function in the IEM response
%         slope2: the slope changing function in the permutation


function [slope1, slope2] = linFitIEM(chanResp, chanRespPermu)

% for IEM channel response 
for f = 1:size(chanResp, 1)

    r = squeeze(chanResp(f,:,:)); %MLW: 6 orientation bins over 188 time points (average over 1D = subjects)
    centerind = 4; %MLW: center of 6 bins is technically 3

    %loop over timepoints and fit the slope ('selectivity') of the iem at each timepoint
    for t = 1:size(r,2) 

        %average the data over equidistant channels
        % why flip the array? -z.x
        avgdat = mean([r(:,t),flipud(r(:,t))],2);
        % why only pick the first half of the channels(from 1 to the center)? -z.x
        avgdat = avgdat(1:centerind);

        %fit a linear model and extract slope and save for each timepoint
        m = polyfit(1:centerind, avgdat', 1);
        slope1(f,t) = m(1); 
    end
end

% for permuation
for f = 1:size(chanRespPermu, 1)

    r = squeeze(chanRespPermu(f,:,:)); %MLW: 6 orientation bins over 188 time points (average over 1D = subjects)
    centerind = 4; %MLW: center of 6 bins is technically 3

    %loop over timepoints and fit the slope ('selectivity') of the iem at each timepoint
    for t = 1:size(r,2) 

        %average the data over equidistant channels
        % why flip the array? -z.x
        avgdat = mean([r(:,t),flipud(r(:,t))],2);
        % why only pick the first half of the channels(from 1 to the center)? -z.x
        avgdat = avgdat(1:centerind);

        %fit a linear model and extract slope and save for each timepoint
        m = polyfit(1:centerind, avgdat', 1);
        slope2(f,t) = m(1); 
    end
end
