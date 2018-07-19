% 2 threshold optimization for IIS detection paper

% load the data frame with all of the sample data
%load(strcat(labDataDrive, '/jonesLab_data/falconHawk/eeg_data/falconHawk_EEGdata_512.mat'));
load('/Users/jesse/Google Drive/IISDetection_Share/Data from Jesse/falconHawk_EEGdata_512.mat');

%%

fs = 512;

% set a master counter to fill the outputoutput matrix
masterNumber = 1;

clear outputoutput

% set threshold levels and find events in the data
highThreshList = [3, 4, 5, 6, 7, 8];
lowThreshList = [-3, -2, -1, 0, 1, 2];

for highthresh = highThreshList
    for lowthresh = lowThreshList
        for i = 1:length(dataFrame)
    
            clear starts ends ideal
            signal = dataFrame(i).signal';
    
            outputoutput.highthesh(masterNumber) = highthresh;
            outputoutput.lowthresh(masterNumber) = lowthresh;
            
            [dataFrame(i).starts, dataFrame(i).ends, dataFrame(i).ideal] = twoThreshPeakDetect(signal, highthresh, lowthresh);
            if length(dataFrame(i).starts) > length(dataFrame(i).ends)
                dataFrame(i).starts = dataFrame(i).starts(1:length(dataFrame(i).ends));
            end

            % merge events that are near to one another
            dataFrame(i).ideal(dataFrame(i).ideal > 0) = 1;

            mergeCuttoff = fs / 5; % 1/5th of a second
            for j = 1:length(dataFrame(i).starts) - 1
                if dataFrame(i).starts(j + 1) - dataFrame(i).ends(j) < mergeCuttoff
                    dataFrame(i).ideal(dataFrame(i).ends(j):dataFrame(i).starts(j+1)) = 1;
                end
            end

            dataFrame(i).starts = [];
            dataFrame(i).ends = [];
            dataFrame(i).starts = find(diff(dataFrame(i).ideal) == 1);
            dataFrame(i).ends = find(diff(dataFrame(i).ideal) == -1);

            if length(dataFrame(i).starts) > length(dataFrame(i).ends)
                dataFrame(i).starts = dataFrame(i).starts(1:end-1);
            end

            % add amplitudes, durations, and build a new ideal signal
            for j = 1:length(dataFrame(i).starts)
                dataFrame(i).eventAmp(j) = max(dataFrame(i).signal(dataFrame(i).starts(j):dataFrame(i).ends(j)));
                dataFrame(i).ideal(dataFrame(i).starts(j):dataFrame(i).ends(j)) = dataFrame(i).eventAmp(j);
                dataFrame(i).durations = dataFrame(i).ends - dataFrame(i).starts;

            end

            dataFrame(i).iei = [0 diff(dataFrame(i).starts)];
            % create a matrix of events that were identified from Matt's two phase 
        end
        clear  eventTreatment eventLabel recordLength eventFilespec

        counter = 1;
        buffer = fs / 5;
        for z = 1:length(dataFrame)
            for j = 1:length(dataFrame(z).starts)
                eventstart = dataFrame(z).starts(j);
                if eventstart > buffer && eventstart < length(dataFrame(z).signal) - buffer
                    eventLabel{counter} = dataFrame(z).uniqueID;
                    eventFilespec{counter} = dataFrame(z).edffilespec;
                    recordLength(counter) = length(dataFrame(z).signal);
                    eventTreatment{counter} = dataFrame(z).treatment;
                    counter = counter + 1;
                end
            end
        end
        
        outputoutput.numEvents(masterNumber) = counter - 1; % counts the number of events identified at each level
        

        summaryTable = table(eventFilespec', eventLabel', eventTreatment', recordLength');
        summaryTable.Properties.VariableNames = {'eventFilespec', 'eventLabel', 'eventTreatment', 'recordLength'};
        summaryTable = grpstats(summaryTable, {'eventTreatment', 'eventLabel', 'eventFilespec', 'recordLength'}, {'sum'});

        summaryTable.GroupCount = summaryTable.GroupCount ./ (summaryTable.recordLength / fs / 3600);
        summaryTable = grpstats(summaryTable, {'eventTreatment', 'eventLabel'}, {'mean'}, 'DataVars', {'GroupCount'});

        summaryTable = summaryTable(~strcmp(summaryTable.eventTreatment, 'WT'), :);
        outputoutput.kaEvents(:,masterNumber) = summaryTable.mean_GroupCount(strcmp(summaryTable.eventTreatment, 'KA'));
        outputoutput.saEvents(:,masterNumber) = summaryTable.mean_GroupCount(strcmp(summaryTable.eventTreatment, 'SA'));
        [~, outputoutput.p(masterNumber),] = ttest2(summaryTable.mean_GroupCount(strcmp(summaryTable.eventTreatment, 'KA')), summaryTable.mean_GroupCount(strcmp(summaryTable.eventTreatment, 'SA')), 'Vartype','unequal');

        disp(strcat({'Round '}, num2str(masterNumber), {' Complete!'}))
        masterNumber = masterNumber + 1;

        
    end
end

%%
save('~/Desktop/twoThresholdOptimization_falconHawk_detectIIS.mat', 'outputoutput')

%%

load('~/Desktop/twoThresholdOptimization_falconHawk_detectIIS.mat')

%%

ticksList = [1, 2, 3, 4, 5, 6];
xtickLabs = {'3', '4', '5', '6', '7', '8'};
ytickLabs = {'-3', '-2', '-1', '0', '1', '2'};
p = outputoutput.p;
highthresh = outputoutput.highthesh;
lowthresh = outputoutput.lowthresh;

figure('units', 'inch', 'pos', [10 10 6 5]) 
pGrid = reshape(p, length(lowThreshList), length(highThreshList));
imagesc(pGrid)
xticks(ticksList)
yticks(ticksList)
xticklabels(xtickLabs)
yticklabels(ytickLabs)
ylabel('low threshold')
xlabel('high threshold')
title('Ttest P value')
colorbar()

caxis([0 0.5])
print('~/Desktop/twoThreshOptimization_ttest.pdf', '-dpdf')

% find min p-value

[minp, minpind] = min(outputoutput.p)

% calculate the rank sum p-value
clear ranksumP
for i = 1:length(outputoutput.p)
    ranksumP(i) = ranksum(outputoutput.kaEvents(:,i), outputoutput.saEvents(:,i));
end


figure('units', 'inch', 'pos', [10 10 6 5]) 
pGrid = reshape(ranksumP, length(lowThreshList), length(highThreshList));
imagesc(pGrid)
xticks(ticksList)
yticks(ticksList)
xticklabels(xtickLabs)
yticklabels(ytickLabs)
ylabel('low threshold')
xlabel('high threshold')
title('ranksum P value')
colorbar()
print('~/Desktop/twoThreshOptimization_ranksumP.pdf', '-dpdf')

%
for i = 1:length(outputoutput.p)
    outputoutput.effectSize(i) = mean(outputoutput.kaEvents(:,i)) - mean(outputoutput.saEvents(:,i));
end

figure('units', 'inch', 'pos', [10 10 6 5]) 
pGrid = reshape(outputoutput.effectSize, length(lowThreshList), length(highThreshList));
imagesc(pGrid)
xticks(ticksList)
yticks(ticksList)
xticklabels(xtickLabs)
yticklabels(ytickLabs)
ylabel('low threshold')
xlabel('high threshold')
title('effectSize')
colorbar()
print('~/Desktop/twoThreshOptimization_effectSize.pdf', '-dpdf')
%

% number of events per animal
for i = 1:length(outputoutput.p)
    
    outputoutput.meanNumEventsPerAnimal(i) = mean([outputoutput.kaEvents(:,i); outputoutput.saEvents(:,i)]);
    
end


figure('units', 'inch', 'pos', [10 10 6 5]) 
pGrid = reshape(outputoutput.meanNumEventsPerAnimal, length(lowThreshList), length(highThreshList));
imagesc(pGrid)
xticks(ticksList)
yticks(ticksList)
xticklabels(xtickLabs)
yticklabels(ytickLabs)
ylabel('low threshold')
xlabel('high threshold')
title('meanNumEventsPerAnimal')
colorbar()
print('~/Desktop/twoThreshOptimization_meanNumEventsPerAnimal.pdf', '-dpdf')

