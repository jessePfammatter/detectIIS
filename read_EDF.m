function EDF = read_EDF(filespec)
% Wrapper for bloackEdfLoad_jp.m that packages everything nicely. This was
% tested to work on xltec EDF files by JP on May 20, 2018 using the KA
% longitudinal data located '/Users/jesse/private/magantiLab_data/CBD_data/EEG files/2017/EF CBD-Torin2 EEG/cohort 4/EDF/Day 1 27 Aug 2017_Animal 2_screenprint.edf'

% Read EDF file
EDF.filespec = filespec;

disp(['Please wait while opening EDF file: ' EDF.filespec])
[EDF.header EDF.signalHeader signalCell] = blockEdfLoad_jp(EDF.filespec);

EDF.signalMat = cell2mat(signalCell);
EDF.signalMax = max( EDF.signalMat );
EDF.signalMin = min( EDF.signalMat );
EDF.signalStd = std( EDF.signalMat );
EDF.lastpt = size(EDF.signalMat, 1);
EDF.nchans = EDF.header.num_signals;
EDF.nrecs = size(EDF.signalMat, 1);
EDF.nsecs = EDF.header.num_data_records;
EDF.fs = EDF.nrecs / EDF.nsecs;

disp('  Done reading EDF')