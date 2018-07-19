% import namefile
namefilespec = strcat(basePath, '/jonesLab_git/jonesLab_code/Projects/falconHawk/falconHawk_namefile_CURRENT.txt');
F = ParseNameFile(namefilespec); % makes the object 'F'

% build a dataFrame file including metadat
n = 1;
for i = 1:length(F.f)
    if str2num(F.itemparams{1, i}.analyzeEEG{1}) == 1
        % add edffilespecs for as many files as listed in the namefile
        for j = 1:length(F.itemparams{1, i}.EEGfile)
            dataFrame(n).edffilespec = strcat(labDataDrive, '/jonesLab_data/falconHawk/eeg_data/', F.itemparams{1, i}.animalGroup{1} , filesep , F.itemparams{1, i}.EEGfile{1, j});
            [a, b, c] = fileparts(dataFrame(n).edffilespec);
            dataFrame(n).uniqueID = strcat(F.itemparams{1, i}.animalGroup{1}, '_', F.itemparams{1, i}.animalID{1});
            dataFrame(n).fs = F.itemparams{1, i}.recordFs{1};
            dataFrame(n).channel = F.itemparams{1, i}.EEGchan{1, j};
            dataFrame(n).treatment = F.itemparams{1, i}.treatment{1};
            n = n + 1;
        end
    end
end

%%
for i = 1:length(dataFrame)

    %actually import and process the data
    alreadyImportedEDF = read_EDF(dataFrame(i).edffilespec);
    dataFrame(i).signal = alreadyImportedEDF.eegData(:,str2num(dataFrame(i).channel));

    % downsample the EEG data to 512 if larger to start with and upsample to 512 if smaller
    if strcmp(dataFrame(i).fs, '1024')
        divisionFactor = 2;
        dataFrame(i).signal = downsample(dataFrame(i).signal, divisionFactor);
        dataFrame(i).fs = str2num(dataFrame(i).fs) / divisionFactor;
    elseif strcmp(dataFrame(i).fs, '256')
        multFactor = 2;
        dataFrame(i).signal = interp(dataFrame(i).signal, multFactor);
        dataFrame(i).fs = str2num(dataFrame(i).fs) * multFactor;
    end
        
    % normalize the signal
    [dataFrame(i).signal, dataFrame(i).model_std, dataFrame(i).modelFit] = normalizeEEG(dataFrame(i).signal, dataFrame(i).fs);
    
end
  

% save the data frame
%save(strcat(labDataDrive, '/jonesLab_data/falconHawk/eeg_data/falconHawk_EEGdata_512.mat'), 'dataFrame', '-v7.3')
save('/Users/jesse/Google Drive/IISDetection_Share/Data from Jesse/falconHawk_EEGdata_512.mat', 'dataFrame', '-v7.3')







