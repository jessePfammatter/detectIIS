function [ output ] = detectEvents_IIS(normSignal, fs)
% DETECTEVENTS_IIS identifies events that contain IISs from a single
% channel EEG record.
%
% This code was originally pulled from 'detectIIS_exploratory_fullData.m'
% on 3/7/18.
%
% JP 2018

    highthresh = 5; % 5 standard deviations above a normalized signal
    lowthresh = -1; % one standard deviation below the mean of a normalized signal.
    
    [starts, ends, ideal] = twoThreshPeakDetect(normSignal', highthresh, lowthresh);
    if size(ends, 2) > 2
        if length(starts) > length(ends)
            starts = starts(1:length(ends));
        end

        % merge events that are near to one another
        ideal(ideal > 0) = 1;

        mergeCuttoff = fs / 5; % half a second
        for j = 1:length(starts) - 1
            if starts(j + 1) - ends(j) < mergeCuttoff
                ideal(ends(j):starts(j+1)) = 1;
            end
        end

        starts = [];
        ends = [];
        starts = find(diff(ideal) == 1);
        ends = find(diff(ideal) == -1);

        if length(starts) > length(ends)
            starts = starts(1:end-1);
        end

        % add amplitudes, durations, and build a new ideal signal
        for j = 1:length(starts)
            eventAmp(j) = max(normSignal(starts(j):ends(j)));
            ideal(starts(j):ends(j)) = eventAmp(j);
            durations = ends - starts;
        end

        iei = [0 diff(starts)];


        % output variables 
        output.starts = starts;
        output.ends = ends;
        output.ideal = ideal;
        output.durations = durations;
        output.amplitudes = eventAmp;
        output.iei = iei;
    else
        output.starts = [];
        output.ends = [];
    end
    

end

