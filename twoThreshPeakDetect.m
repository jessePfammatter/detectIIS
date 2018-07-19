function [starts, ends, ideal] = twoThreshPeakDetect(signal, highthresh, lowthresh)
% Matt's two threshold peak detect
%
% modified by JP 2017 

% 

abovehigh = 0.*signal;
abovehigh(find(signal>highthresh)) = 1; % changed to 1 from highthresh by JP to allow values at 0 to work.
if sum(abovehigh) > 0
    highcross = [0 diff(abovehigh)];        % prepend with a zero to preserve the number of points and alignment when using diff()
    highcross(highcross<0) = 0;             % delete any negative crossings of this threshold

    % Find all points below lowthresh and mark only thresh crossings in the negative direction
    belowlow = 0.*signal;
    belowlow(find(signal<lowthresh)) = 1;   % changed to 1 from highthresh by JP to allow values at 0 to work.
    lowcross = [0 diff(belowlow)];          % prepend with a zero to preserve the number of points and alignment when using diff()
    lowcross(lowcross<0) = 0;               % delete any positive crossings of this threshold

    % To ensure that the first detected event is a highcrossing, find any lowcrosses
    % that occur before the first highcross and delete them
    dummy1 = find(highcross);
    lowcross(1:dummy1(1)) = 0;

    % Combine highcrosses and lowcrosses to give an "idealized" signal that is
    % zero until a highcross (1) and stays that way until a lowcross (0)
    % This is probably the slowest part. There must be a better way than using loops and ifs. Vectorizing this
    % somehow might be important when dealing with millions of data points in
    % 24 hour EEG records or whatever
    ideal = 0.*signal;
    state = 0;
    for n = 1:length(signal)
        if highcross(n)
            state = 1;
        elseif lowcross(n)
            state = 0;
        end
        ideal(n) = highthresh.*state;
    end    

    % Finally, find the actual start and end times (points, in this example) of each detected "event"
    dummy = diff( sign([0 ideal]));     % prepend with a zero to preserve the number of points and alignment when using diff()
    starts  = find(dummy==1);
    ends    = find(dummy==-1);
else
    starts = [];
    ends = [];
    ideal = abovehigh;
end

end