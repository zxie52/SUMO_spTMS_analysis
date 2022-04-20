%% Quantify the selectivity of IEM reconstruction via linear fit

% Notes:
% I made up data here, but r should be the centered reconstruction averaged 
% across trials for a single subject. You can apply the same code to 
% single-trial reconstructions but it will probably be very noisy.
%
% If the number of orientation channels (here 9) is even, then you'll need
% to copy the response from the channel furthest from the center
% and paste it on the other end of the matrix to make the first dim
% of r odd. See example at the bottom of this script for how to do this.
%
% The center of r (e.g., r(5,:) needs to be the center of the tuning curve. 

types = {chanresps1LA, chanresps1LU, chanresps1RA, chanresps1RU};
subjects = {100, 102, 103, 104, 105, 106, 108, 109};    

%chanresps2LA_stay, chanresps2LA_switch, chanresps2LU_stay, chanresps2LU_switch,chanresps2RA_stay, chanresps2RA_switch, chanresps2RU_stay, chanresps2RU_switch};

for f = 1:length(subjects)
%for e = 1:length(types)

    %made up data: 9-channel reconstructions over 200 time points
    %r = randn(9,200);
    r = squeeze(chanresps1RU(f,:,:)); %MLW: 6 orientation bins over 188 time points (average over 1D = subjects)
    centerind = 3; %MLW: center of 6 bins is technically 3

    %specifiy number of orientation information channels in the iem, here 9, must be odd
    %nchan = 9;
    nchan = 7; %MLW20180716: we have 6 orientation bins, but 7 after we duplicate

    % if you have an even number of channels in your reconstruction you can
    % make it odd by duplicating the channel furthest off-center prior to
    % fitting a slope - i.e we simply do the following before fitting a slope via the code above
    % find furthest channel and dupliate it in the appropriate position
    if length(centerind:size(r,1)) > centerind
        chan2dup = size(r,1);
        r = [r(chan2dup,:);r];
    else
        chan2dup = 1;
        r = [r;r(chan2dup,:)];
    end

    % now the new center channel is the middle of r, for example.
    newcent = ceil(size(r,1)/2);

    %loop over timepoints and fit the slope ('selectivity') of the iem at each timepoint
    for t = 1:size(r,2) 

        %average the data over equidistant channels
        % why flip the array? -z.x
        avgdat = mean([r(:,t),flipud(r(:,t))],2);
        % why only pick the first half of the channels(from 1 to the center)? -z.x
        avgdat = avgdat(1:ceil(nchan/2));

        %fit a linear model and extract slope and save for each timepoint
        m = polyfit(1:ceil(nchan/2), avgdat', 1);
        slope(f,t) = m(1); 

    end

%end
end
    
% Now we can plot our slope estimates over time. This is made up data so
% it should just fluctuate around zero. Larger positive values mean greater
% selectivity, negative values mean negetive selectivity (whatever that means)
figure; plot(dstimes, squeeze(mean(slope,1)));
xlabel('time from stimulus onset');
ylabel('slope');
xlim([dstimes(1) dstimes(end)]);

% Because this is data from just one subject and we fit the
% trial-averaged recontruction, we have no estimate of variability.
% You can repeat the above steps for each subject, and then plot the
% average slope over time with a standard error or 95% CI. you can test
% for significane by comparing the slope against zero at each timepoint
% with a t-test, or by using non-parametric bootstrapping (as in our JoCN paper)
% you can use bootstrapping to assess single-subject significance as well
% let me know if you wanna go that route

%% bonus material

% if you have an even number of channels in your reconstruction you can
% make it odd by duplicating the channel furthest off-center prior to
% fitting a slope.

% for example if my reconstruction is 8 channels, and the 4th index is the  center channel
r = randn(8,200);
centerind = 4;

% then we can simply do the following before fitting a slope via the code above.
% find furthest channel and dupliate it in the appropriate position
if length(centerind:size(r,1)) > centerind
    chan2dup = size(r,1);
    r = [r(chan2dup,:);r];
else
    chan2dup = 1;
    r = [r;r(chan2dup,:)];
end

% now the new center channel is the middle of r, for example.
newcent = ceil(size(r,1)/2);

%% extra special bonus material

% Let's convince ourselves that the slope fitting procedure worked. Here
% we'll introduce a robust tuning curve at from 110 to 140 samples into the trial.

%random data
r = randn(9,200);

%create cosine signal
sig = cos(linspace(-pi,pi,size(r,1)))*2; %2 here controls the magnitude of the cos;

%add signal to part of the data;
r(:,110:140) = r(:,110:140) + repmat(sig,140-110+1,1)';

%now use this r new in the slope fitting code above;
