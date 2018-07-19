function [ normSignal, sig, modelfit, mu] = normalizeEEG( signal , fs)
% normalizeEEG normalizes an EEG signal using the mu and sigma parameters
% from a guassian curve fit the the all points histogram of the signal. The
% resultant normalized signal is similar to a Z-score normalization but
% produces better results in the cases of singal noise or significant
% amounts of high amplitdue signal (e.g. seizures).
%
% [ normSignal, sig, modelfit, mu] outputs provide:
%
% normSignal: normalized signal
% sig: the standard deviation of the model fit
% modelFit: a output describing how well the guassian fit the data where 1
% is a perfect fit and 0 is the worst fit on the planet yourDataSucks.
% mu: the mean of the guassian fit
%
% example:
%
% %edffilespec = '/Users/jesse/private/CBD_Torin2/EEG files/2017/Kcna_CBD/Set_1/Day 1/KCNA_CBD_Dec 2017_Cohort 1_Day 1_An1.edf';
% edffilespec = '/Volumes/cookieMonster/Kcna_CBD/Set_1/Day 1/KCNA_CBD_Dec 2017_Cohort 1_Day 1_An1.edf';
% [header] = edfread(edffilespec);
% targetSignal = find(strcmp(header.label, strcat('RF', num2str( 1)))); % looking for channel RF1 for this animal
% [header, signal] = edfread(edffilespec,'targetSignals',targetSignal);
% fs = header.frequency(1);
% [normSignal, sig, modelfit, mu] = normalizeEEG( signal , fs);
%
% JP 2017, modified 2018

% apply filter to signal
filtSignal = sixtyHzFilt_EEG(signal, fs); % notch filter at 60Hz
filtSignal = highPassChebyshev1Filt_EEG(filtSignal, fs); % high pass filter above 0.5Hz

% remove stretches of flatline zeros in the record. This does not influence the analysis if it takes out a few actual zeros (there are few of them in a quality recording.
tempSignal = filtSignal(filtSignal ~= 0);

% fit a guassian model to the signal
hfit = histfit(tempSignal, floor(sqrt(length(tempSignal))), 'normal');
y = hfit(1).YData;
x = hfit(1).XData;
guessInd = find(y == max(y));
guess = [x(guessInd), std(tempSignal), max(y)];
[guess, fval] = fminsearch( 'fit_gauss', guess, [], x, y, 1);
hold on;
[mu, sig, amp] = deal(  guess(1), guess(2), guess(3) ); 
est = amp.*exp(-(x-mu).^2./sig.^2);
bar(x, y);
plot(x, est, 'g');
title('Hist/Fit', 'interpreter', 'none');
ylabel('sqrt(length(EEGsignal))');
xlabel('Current (pA)');
hold off;

% calculate model fit. 
auc = trapz(x, est); % area under model
aud = trapz(x, y); % area under data (usually the larger number)

modelfit = 1- ((aud - auc) / aud); % values less than one means that data area is larger than the model fit. This is normal since the model should fit a portion of the data and extreme values should be 'outside' the model fit.

% adjust the signal by the mean and sigma value of the model fit.
tempSignal = filtSignal - mu;
tempSignal = tempSignal / sig;
normSignal = tempSignal;

% close window after plot
close all

end

