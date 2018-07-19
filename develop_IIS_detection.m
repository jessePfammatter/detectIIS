%% STARTUP FOR IIS DETECTION. RUN THIS FIRST SECTION AND THEN EITHER ONE OF THE SECTIONS BELOW

% load the falconHawk data
%load(strcat(labDataDrive, '/jonesLab_data/falconHawk/eeg_data/currentWorkingFile.mat'));
load('/Users/jesse/Google Drive/research/jonesLab/manuscripts/IIS detection manuscript/IISDetection_Share/Data from Jesse/falconHawk_EEGdata_512.mat');
disp('All Done!')

%% skip if happy with all of the variables in the load file above then don't need to run this section

% check the sample rate
for i = 1:length(dataFrame)
    fsList(i) = dataFrame(i).fs;
end

if length(unique(fsList)) == 1
    fs = unique(fsList);
else
    error('Different sample rates for each file. Adjust accordingly.')
end


% remove the event amplitude field if it exists
if isfield(dataFrame, 'eventAmp')
    dataFrame = rmfield(dataFrame, 'eventAmp');
end
    
for i = 1:length(dataFrame)
    
    clear starts ends ideal
    signal = dataFrame(i).signal';
    
    % set threshold levels and find events in the data, this was picked as the lowest p-value from the threshold optimization data set.
    highthresh = 5;
    lowthresh = -1;
    [dataFrame(i).starts, dataFrame(i).ends, dataFrame(i).ideal] = twoThreshPeakDetect(signal, highthresh, lowthresh);
    if length(dataFrame(i).starts) > length(dataFrame(i).ends)
        dataFrame(i).starts = dataFrame(i).starts(1:length(dataFrame(i).ends));
    end
 
    % merge events that are near to one another
    dataFrame(i).ideal(dataFrame(i).ideal > 0) = 1;

    mergeCuttoff = fs / 5; % 200 ms
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
end

%
% create a matrix of events that were identified from Matt's two phase 
clear eventMatrix eventStart eventTreatment eventLabel eventAmplitude eventDuration recordLength eventFilespec

for i = 1:length(dataFrame)
    for j = 1:length(dataFrame(i).starts)
    end
end
counter = 1;
buffer = fs; % perhaps switch to fs for 1 second on each side?
for i = 1:length(dataFrame)
    for j = 1:length(dataFrame(i).starts)
        eventstart = dataFrame(i).starts(j);
        if eventstart > buffer && eventstart < length(dataFrame(i).signal) - buffer
            eventMatrix(counter, :) = dataFrame(i).signal(eventstart-buffer:eventstart+buffer);
            eventStart(counter, :) = eventstart;
            eventLabel{counter} = dataFrame(i).uniqueID;
            eventFilespec{counter} = dataFrame(i).edffilespec;
            recordLength(counter) = length(dataFrame(i).signal);
            eventTreatment{counter} = dataFrame(i).treatment;
            eventDuration(counter) = dataFrame(i).durations(j);
            eventAmplitude(counter) = dataFrame(i).eventAmp(j);
            counter = counter + 1;
        end
    end
end

%{
built this for the KcnA analysis stuff. but also this is a good file to
have..

KA_out.eventMatrix = eventMatrix;
KA_out.eventDurations = eventDuration;
KA_out.eventStart = eventStart;
KA_out.eventFilespec = eventFilespec;
KA_out.animalID = eventLabel;
KA_out.recordLength = recordLength;
KA_out.eventTreatment = eventTreatment;

% save this eventMatrix 
save('/Volumes/cookieMonster/Kcna_CBD/KA_eventMatrix.mat', 'KA_out', '-v7.3') 


% make a color list for kainate vs saline
for i = 1:length(dataFrame)
    if strcmp(dataFrame(i).treatment, {'KA'})
        colorcolor{i} = 'r';
    else
        colorcolor{i} = 'k';
    end
end
%}

% EPILEPSY INDEX FROM CLUSTERING METHOD prep stuff

% some prep variable
options = statset('MaxIter',100000);

% total number of saline and kainate recoriding hours in the file.
totalDataPointsSA = 0;
totalDataPointsKA = 0;
for i = 1:length(dataFrame)
    if strcmp(dataFrame(i).treatment, {'SA'})
        totalDataPointsSA = totalDataPointsSA + length(dataFrame(i).signal);
    elseif strcmp(dataFrame(i).treatment, {'KA'})
        totalDataPointsKA = totalDataPointsKA + length(dataFrame(i).signal);
    end
end
totalHoursSA = totalDataPointsSA / fs / 3600;
totalHoursKA = totalDataPointsKA / fs / 3600;

% cluster optimization
clear optimizationOutput
counter = 1;
for randomSeed = 1:25
    for nClusters = 2:25
        for nPCAComponents = [3] % how many pca components matter? % try 3, 6, 9
            for dataInputMatrix = 2 % options are 1:2
                for addDurationData = 1 % options are 1:2
                %{
                    % PCA for variable selection
                    if dataInputMatrix == 1
                        PCA = pcamj(waveMat , 1); % viewValues = [-113 15]; xAxisValues = [-200 100]; yAxisValues = [-100 100]; zAxisValues = [-80 60];
                    else
                        PCA = pcamj(eventMatrix , 1); %viewValues = [-110 35]; xAxisValues = [-150 150]; yAxisValues = [-100 100]; zAxisValues = [-80 60];
                    end
%}
                    
                    PCA = pcamj(eventMatrix , 1); %viewValues = [-110 35]; xAxisValues = [-150 150]; yAxisValues = [-100 100]; zAxisValues = [-80 60];
                                            
                    % what are the optimization settings
                    optimizationOutput(counter).randomSeed = randomSeed;
                    optimizationOutput(counter).nClusters = nClusters;
                    optimizationOutput(counter).nPCAComponents = nPCAComponents;
                    optimizationOutput(counter).dataInputMatrix = dataInputMatrix;
                    optimizationOutput(counter).addDurationData = addDurationData;

                    % run the gm and clustering
                    rng(randomSeed); % For reproducibility and random seed checking
                    pcadata = PCA.Proj(:,1:nPCAComponents);
                    
                    if addDurationData == 2
                        pcadata = [pcadata, eventDuration'];
                    end
                    
                    optimizationOutput(counter).gm = fitgmdist(pcadata,nClusters,...
                                                'CovarianceType','full',...
                                                'SharedCovariance',false,...
                                                'RegularizationValue',.01,...
                                                'Options',options,...
                                                'Start', 'plus');
                                            
                    [optimizationOutput(counter).clusterX optimizationOutput(counter).NLOGL] = cluster(optimizationOutput(counter).gm, pcadata);
                    P = posterior(optimizationOutput(counter).gm,pcadata);
                    temptemp = sum((P.^2)');
                    minVal = 1/nClusters;
                    temptemp = temptemp .* (1 - ((1 - temptemp) .* minVal));
                    optimizationOutput(counter).clusterConfidenceMeas_mean = mean(temptemp);
                    optimizationOutput(counter).clusterConfidenceMeas_std = std(temptemp);

                    % increment counter
                    counter = counter + 1;
                    counter - 1 % to diplsay which iteration we are on
                end
            end
        end
    end
end

% calculate the cluster probabilities
for counter = 1:length(optimizationOutput)
    idx = optimizationOutput(counter).clusterX;
    optimizationOutput(counter).clusterEpilepsyRatio_odds = 0;
    optimizationOutput(counter).clusterEpilepsyRatio_prob = 0;

    for j = 1:2
        if j == 1
            treatmentGroup = 'SA';
            treatmentColor = {'k'};
        else
            treatmentGroup = 'KA';
            treatmentColor = {'r'};
        end
        
        for i = 1:max(idx)
            temp = eventMatrix(intersect(find(idx == i), find(strcmp(eventTreatment, {treatmentGroup}))),fs-128:fs+128); 
            numElements = size(temp, 1);
            if j == 1
                clusterEpilepsyRatioSA(i) = numElements / totalHoursSA;
            else
                clusterEpilepsyRatioKA(i) = numElements / totalHoursKA;
            end
        end
    end

    optimizationOutput(counter).clusterEpilepsyRatio_odds = clusterEpilepsyRatioKA ./ clusterEpilepsyRatioSA; % the odds ratio of being an event related to KA.
    optimizationOutput(counter).clusterEpilepsyRatio_prob = optimizationOutput(counter).clusterEpilepsyRatio_odds ./  (optimizationOutput(counter).clusterEpilepsyRatio_odds + 1); % 
    optimizationOutput(counter).clusterEpilepsyRatio_prob(find(isnan(optimizationOutput(counter).clusterEpilepsyRatio_prob))) = 1; % 
    clear clusterEpilepsyRatioKA clusterEpilepsyRatioSA

end

% make a list of all of the uniqueIDs and the total number of hours of recording for each of those
for i = 1:length(dataFrame)
    recordLengths(i) = length(dataFrame(i).signal) / fs / 3600;
    animalName{i} = dataFrame(i).uniqueID;
end

recordLengthsTable = table(recordLengths', animalName');
recordLengthsTable.Properties.VariableNames = {'recordLengths', 'animalName'};
recordLengthsTable = grpstats(recordLengthsTable, {'animalName'}, {'sum'}, 'DataVars',{'recordLengths'});
    
% calculate the epiepsy index
clear probsList
for counter = 1:length(optimizationOutput)
    idx = optimizationOutput(counter).clusterX;
    for i = 1:max(idx)
        probsList(idx == i) = optimizationOutput(counter).clusterEpilepsyRatio_prob(i);
    end
    summaryTable = table(idx, probsList', eventLabel', eventTreatment');
    summaryTable.Properties.VariableNames = {'idx', 'probsList', 'eventLabel', 'eventTreatment'};
    summaryTable = grpstats(summaryTable, {'eventTreatment', 'eventLabel', 'idx'}, {'mean'}, 'DataVars',{'probsList'});
    summaryTable.eventsPerHour = summaryTable.GroupCount;
    for i = 1:size(summaryTable, 1)
        summaryTable.eventsPerHour(i) = summaryTable.GroupCount(i) / recordLengthsTable.sum_recordLengths(strcmp(summaryTable.eventLabel(i), recordLengthsTable.animalName));
    end

    summaryTable.mean_probsList = (summaryTable.mean_probsList - .5) .* 2;
    summaryTable = summaryTable(summaryTable.mean_probsList > 0,:);
    summaryTable.epilepsyEventsPerHour = summaryTable.eventsPerHour .* summaryTable.mean_probsList;

    summaryTable = grpstats(summaryTable, {'eventTreatment', 'eventLabel'}, {'sum'}, 'DataVars',{'epilepsyEventsPerHour'});
    summaryTable = summaryTable(~strcmp(summaryTable.eventTreatment, 'WT'), :);

    [~, optimizationOutput(counter).p,] = ttest2(summaryTable.sum_epilepsyEventsPerHour(strcmp(summaryTable.eventTreatment, 'KA')), summaryTable.sum_epilepsyEventsPerHour(strcmp(summaryTable.eventTreatment, 'SA')), 'Vartype','unequal');
    summaryTable = grpstats(summaryTable, {'eventTreatment'}, {'mean', 'meanci'}, 'DataVars',{'sum_epilepsyEventsPerHour'});
    optimizationOutput(counter).kamean = summaryTable.mean_sum_epilepsyEventsPerHour(strcmp(summaryTable.eventTreatment, 'KA'));
    optimizationOutput(counter).samean = summaryTable.mean_sum_epilepsyEventsPerHour(strcmp(summaryTable.eventTreatment, 'SA'));
    optimizationOutput(counter).kaci = summaryTable.meanci_sum_epilepsyEventsPerHour(strcmp(summaryTable.eventTreatment, 'KA'),:);
    optimizationOutput(counter).saci = summaryTable.meanci_sum_epilepsyEventsPerHour(strcmp(summaryTable.eventTreatment, 'SA'),:);
    
end

disp('All Done!')


%% Plots that before selection of high and low threshold -- optimization stuff and events per hour at whichever two thresholds we've selected
% now plot the optimization parameters
clear dataInputMatrix nClusters nPCAComponents effectSize p saci clusterConfidenceMeas_mean clusterConfidenceMeas_std NLOGL addDurationData

for counter = 1:length(optimizationOutput)
    
    dataInputMatrix(counter) = optimizationOutput(counter).dataInputMatrix; % wavelet or raw data
    nClusters(counter) = optimizationOutput(counter).nClusters; % the number of clusters
    nPCAComponents(counter) = optimizationOutput(counter).nPCAComponents;
    effectSize(counter) = optimizationOutput(counter).kamean - optimizationOutput(counter).samean;
    p(counter) = optimizationOutput(counter).p;
    saci(counter) = optimizationOutput(counter).saci(2);
    addDurationData(counter) = optimizationOutput(counter).addDurationData;
    clusterConfidenceMeas_mean(counter) = optimizationOutput(counter).clusterConfidenceMeas_mean;
    if clusterConfidenceMeas_mean(counter) > 1
        clusterConfidenceMeas_mean(counter) = 1;
    end
    clusterConfidenceMeas_std(counter) = optimizationOutput(counter).clusterConfidenceMeas_std;
    NLOGL(counter) = optimizationOutput(counter).NLOGL;
end

% ----- for Figure 1 events per hour index.
h = figure('units', 'inch', 'pos', [10 10 5 5]) 

% create summary table for the figure
summaryTable = table(eventLabel', eventTreatment', recordLength', eventFilespec');
summaryTable.Properties.VariableNames = {'eventLabel', 'eventTreatment', 'recordLength', 'eventFilespec'};
summaryTable = grpstats(summaryTable, {'eventTreatment', 'eventLabel', 'recordLength', 'eventFilespec'}, {'sum'});
summaryTable.GroupCount = summaryTable.GroupCount ./ (summaryTable.recordLength / fs / 3600);
summaryTable = grpstats(summaryTable, {'eventTreatment', 'eventLabel'}, {'mean'}, 'DataVars', {'GroupCount'});
summaryTable = summaryTable(~strcmp(summaryTable.eventTreatment, 'WT'), :);
[~, ppp,] = ttest2(summaryTable.mean_GroupCount(strcmp(summaryTable.eventTreatment, 'KA')), summaryTable.mean_GroupCount(strcmp(summaryTable.eventTreatment, 'SA')), 'Vartype','unequal')
boxplot(summaryTable.mean_GroupCount, summaryTable.eventTreatment, 'color', 'k')
hold on
for i = 1:size(summaryTable.eventTreatment, 1)
    if strcmp(summaryTable.eventTreatment(i), 'KA')
        if strcmp(summaryTable.eventLabel(i), '1-5b_L1R0') || strcmp(summaryTable.eventLabel(i), '1-6b_L1R2')
            plot(normrnd(1, .05), summaryTable.mean_GroupCount(i), 'm.', 'MarkerSize', 25) %
        else
            plot(normrnd(1, .05), summaryTable.mean_GroupCount(i), 'ro', 'MarkerSize', 7)
        end
    else
        plot(normrnd(2, .05), summaryTable.mean_GroupCount(i), 'ko', 'MarkerSize', 7)
    end
end
ylabel('Events Per Hour')
xlabel('Treatment')
text(1.5, 80, num2str(ppp))

% paper settings and print
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
print('~/Desktop/eventsPerHour_twoThresh.pdf', '-dpdf')

% ----- Plots and summary stuff for optimization (Figure 3)
h = figure('units', 'inch', 'pos', [10 10 12 7]) 
subplot(2,2,1)
    scatter(nClusters, NLOGL, 10, 'k', 'filled')
    xlabel('Number of Clusters')
    ylabel('Negative Loglikelihood')
    breakyaxis([3.6e5, 4.35e5]) 
subplot(2,2,2)
    scatter(nClusters, clusterConfidenceMeas_mean, 10, 'k', 'filled')
    xlabel('Number of Clusters')
    ylabel('Cluster Confidence Measure')
subplot(2,2,3)
    scatter(nClusters, effectSize, 10, 'k', 'filled')
    xlabel('Number of Clusters')
    ylabel('Effect Size of Epilepsy Index')
subplot(2,2,4)
    scatter(nClusters, p, 10, 'k', 'filled')
    xlabel('Number of Clusters')
    ylabel('P-value for Epilepsy Index')

% paper settings and print
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])    
print('~/Desktop/figure3_detectIIS.pdf', '-dpdf')


% plot the variance explained by the PCA, using this for anything ate % variance explained at 3 PC
PCA = pcamj(eventMatrix, 1);
varianceComps = diag(fliplr(PCA.D));
varianceComps = varianceComps ./ max(varianceComps);
perVarComps = flipud((1 - varianceComps) * 100); 

h = figure('units', 'inch', 'pos', [10 10 5 5]) 
    plot(perVarComps)
    title('rawMatrix')
    ylim([0, 100])
    xlim([0, 150])
    ylabel('% of Variance Explained by PCs')
    xlabel('Principle Componenets')
    text(50, 50, num2str(perVarComps(3)))
    
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
print('~/Desktop/PCAVarianceComponents_detectIIS.pdf', '-dpdf')

% for figure 1 A, B
windowwindow = (1:10000000) + 15505000;
windowwindowTime = (windowwindow / fs / 3600) + 8;
h = figure('units', 'inch', 'pos', [10 10 25 10]);
subplot(3, 2, 1)
    whichFile = 19;
    hold on;
    plot(windowwindowTime , dataFrame(whichFile).signal(windowwindow), 'k')
    plot(windowwindowTime , dataFrame(whichFile).ideal(windowwindow), 'r')
    xlim([windowwindowTime(1) , windowwindowTime(end) ]);
    ylim([-10 25])
    xlabel('Time of Day (h)')
    title('KA Animal')
subplot(3, 2, 5)
    whichFile = 19;
    zoomWindow = windowwindow(960000:965000) ;
    zoomWindowTime = (zoomWindow / fs / 3600) + 8;
    hold on;
    plot(zoomWindowTime, dataFrame(whichFile).signal(zoomWindow), 'k') 
    plot(zoomWindowTime, dataFrame(whichFile).ideal(zoomWindow), 'r')
    xlabel('Time of Day (h)')
    xlim([zoomWindowTime(1) , zoomWindowTime(end) ]);
    ylim([-10 20])   
subplot(3, 2, 3)
    whichFile = 19;
    bumperbumper = 900000;
    zoomWindow = windowwindow(960000-bumperbumper:965000+bumperbumper) ;
    zoomWindowTime = (zoomWindow / fs / 3600) + 8;
    hold on;
    plot(zoomWindowTime, dataFrame(whichFile).signal(zoomWindow), 'k') 
    plot(zoomWindowTime, dataFrame(whichFile).ideal(zoomWindow), 'r')
    xlabel('Time of Day (h)')
    xlim([zoomWindowTime(1) , zoomWindowTime(end) ]);
    ylim([-10 20])       
subplot(3, 2, 1)
    plot(zoomWindowTime(1), -8, '.b')   
subplot(3, 2, 2)
    whichFile = 5;
    hold on;
    plot(windowwindowTime, dataFrame(whichFile).signal(windowwindow), 'k')
    plot(windowwindowTime , dataFrame(whichFile).ideal(windowwindow), 'r')
    ylim([-10 25])
    xlabel('Time of Day (h)')
    xlim([windowwindowTime(1) , windowwindowTime(end) ]);
    title('SA Animal')
subplot(3, 2, 6)
    whichFile = 5;
    zoomWindow = windowwindow(2010500:2015500);
    zoomWindowTime = (zoomWindow / fs / 3600) + 8;
    hold on
    plot(zoomWindowTime , dataFrame(whichFile).signal(zoomWindow), 'k')
    plot(zoomWindowTime, dataFrame(whichFile).ideal(zoomWindow), 'r')
    xlabel('Time of Day (h)')
    xlim([zoomWindowTime(1) , zoomWindowTime(end) ]);
    ylim([-10 20])    
subplot(3, 2, 4)
    whichFile = 5;
    zoomWindow = windowwindow(2010500-bumperbumper:2015500+bumperbumper);
    zoomWindowTime = (zoomWindow / fs / 3600) + 8;
    hold on
    plot(zoomWindowTime , dataFrame(whichFile).signal(zoomWindow), 'k')
    plot(zoomWindowTime, dataFrame(whichFile).ideal(zoomWindow), 'r')
    xlabel('Time of Day (h)')
    xlim([zoomWindowTime(1) , zoomWindowTime(end) ]);
    ylim([-10 20])


set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
print('~/Desktop/figure1_IIS.pdf', '-dpdf', '-painters')

close all

%%
% find the iterations with 3 PC, 9 clusters, raw matrix with the lowest log likelihood -- the optimial model to make the next set of plots with.

clear dataInputMatrix nClusters nPCAComponents effectSize p saci clusterConfidenceMeas_mean clusterConfidenceMeas_std NLOGL

for counter = 1:length(optimizationOutput)
    
    dataInputMatrix(counter) = optimizationOutput(counter).dataInputMatrix; % wavelet or raw data
    nClusters(counter) = optimizationOutput(counter).nClusters; % the number of clusters
    nPCAComponents(counter) = optimizationOutput(counter).nPCAComponents;
    effectSize(counter) = optimizationOutput(counter).kamean - optimizationOutput(counter).samean;
    p(counter) = optimizationOutput(counter).p;
    saci(counter) = optimizationOutput(counter).saci(2);
    clusterConfidenceMeas_mean(counter) = optimizationOutput(counter).clusterConfidenceMeas_mean;
    if clusterConfidenceMeas_mean(counter) > 1
        clusterConfidenceMeas_mean(counter) = 1;
    end
    clusterConfidenceMeas_std(counter) = optimizationOutput(counter).clusterConfidenceMeas_std;
    NLOGL(counter) = optimizationOutput(counter).NLOGL;
end

index = nPCAComponents == 3 & nClusters == 9 & dataInputMatrix == 2;
[a, b] = min(NLOGL(index));
useThisOne = find(NLOGL == a);

% rebuild important objects from optimized event
rng(optimizationOutput(useThisOne).randomSeed); % For reproducibility
nClusters = optimizationOutput(useThisOne).nClusters;
nPCAComponents = optimizationOutput(useThisOne).nPCAComponents;
gmBest = optimizationOutput(useThisOne).gm;
idx = optimizationOutput(useThisOne).clusterX;
PCA = pcamj(eventMatrix, 1);
pcadata = PCA.Proj(:,1:nPCAComponents);

% calculate posterior probabilities
P = posterior(gmBest,pcadata); 
n = size(pcadata,1);
[~,order] = sort(P(:,2));


%%
% plot the pca (Figure 2)
% this is the 3d plot similar to the first version of the mansucript. Below is a set of 4 2d plots. 
h = figure('units', 'inch', 'pos', [25 10 17 7]) 

subplot(1,2,1)
    scatter3(pcadata(strcmp(eventTreatment, 'KA'),1),pcadata(strcmp(eventTreatment, 'KA'),2), pcadata(strcmp(eventTreatment, 'KA'),3), 10, 'r', 'filled')
    hold on;
    scatter3(pcadata(strcmp(eventTreatment, 'SA'),1),pcadata(strcmp(eventTreatment, 'SA'),2), pcadata(strcmp(eventTreatment, 'SA'),3), 10, 'k', 'filled')
    grid on;
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    xlim([-275 200])
    ylim([-100 150])
    zlim([-120 100])
    view(67,-32) 

subplot(1,2,2)
    scatter3(pcadata(:,1), pcadata(:,2), pcadata(:,3), 10, idx, 'filled')
    colormap('redgreenblue')
    grid on;
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')

    xlim([-275 200])
    ylim([-100 150])
    zlim([-120 100])
    view(67,-32)
 
   %{
    % for alternate view
    xlim([-350 250])
    ylim([-250 200])
    zlim([-180 160])
    view(-106, 56) 
    %}

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
print(strcat('~/Desktop/clusterSummaryA1.pdf'), '-dpdf', '-painters')
%print(strcat('~/Desktop/clusterSummaryA1alt.pdf'), '-dpdf', '-painters')
close all

%{
% sets of 2d plots that describe the space

figure('units', 'inch', 'pos', [25 10 12 15]) 
markerAlphaLevel = 0.2;
subplot(2,2,1)
    scatter(pcadata(strcmp(eventTreatment, 'KA'),1),pcadata(strcmp(eventTreatment, 'KA'),2), 10, 'r', 'filled')
    hold on;
    scatter(pcadata(strcmp(eventTreatment, 'SA'),1),pcadata(strcmp(eventTreatment, 'SA'),2), 10, 'k', 'filled')
    grid on;
    xlabel('PC1')
    ylabel('PC2')
    xlim([-200 300])
    ylim([-100 150])

subplot(2,2,2)
    scatter(pcadata(strcmp(eventTreatment, 'KA'),3),pcadata(strcmp(eventTreatment, 'KA'),2), 10, 'r', 'filled')
    hold on;
    scatter(pcadata(strcmp(eventTreatment, 'SA'),3),pcadata(strcmp(eventTreatment, 'SA'),2), 10, 'k', 'filled')
    grid on;
    xlabel('PC3')
    ylabel('PC2')
    xlim([-120 120])
    ylim([-100 150])
    
subplot(2,2,3)
    scatter(pcadata(:,1), pcadata(:,2), 10, idx, 'filled')
    colormap('redgreenblue')
    grid on;
    xlabel('PC1')
    ylabel('PC2')
   xlim([-200 300])
    ylim([-100 150])

subplot(2,2,4)
    scatter(pcadata(:,3), pcadata(:,2), 10, idx, 'filled')
    colormap('redgreenblue')
    grid on;
    xlabel('PC3')
    ylabel('PC2')
    xlim([-120 120])
    ylim([-100 150])

    print(strcat('~/Desktop/clusterSummaryA2.pdf'), '-dpdf', '-painters')

%}

% plot the ensemble of events from each cluster and a bunch of other things from the clustering
totalDataPointsSA = 0;
totalDataPointsKA = 0;
for i = 1:length(dataFrame)
    if strcmp(dataFrame(i).treatment, {'SA'})
        totalDataPointsSA = totalDataPointsSA + length(dataFrame(i).signal);
    elseif strcmp(dataFrame(i).treatment, {'KA'})
        totalDataPointsKA = totalDataPointsKA + length(dataFrame(i).signal);
    end
end
totalHoursSA = totalDataPointsSA / fs /3600;
totalHoursKA = totalDataPointsKA / fs /3600;
clear clusterEpilepsyRatio_odds clusterEpilepsyRatio_prob clusterEpilepsyRatioKA clusterEpilepsyRatioSA

%plot the number of events in each cluster that are KA vs. Saline that belong to each cluster 
for j = 1:2
    if j == 1
        treatmentGroup = 'SA';
        treatmentColor = {'k'};
    else
        treatmentGroup = 'KA';
        treatmentColor = {'r'};
    end
   h = figure('units', 'inch', 'pos', [10 10 15 15]) 
    for i = 1:max(idx)
        subplot(3, 3, i)
            temp = eventMatrix(intersect(find(idx == i), find(strcmp(eventTreatment, {treatmentGroup}))),fs-128:fs+128); 
            numElements = size(temp, 1);
            plot(temp')
            patchVectorA = mean(temp) + (std(temp));
            patchVectorB = mean(temp) - (std(temp));
            hold on;
            lengthLength = size(patchVectorA, 2);
            if j == 1
                patch([1:lengthLength, fliplr(1:lengthLength)], [patchVectorA, fliplr(patchVectorB)], 'k', 'FaceColor', 'k', 'EdgeColor', 'none', 'facealpha', .4) % , colorcolor, 'FaceColor', colorcolor, 'EdgeColor', 'none'
                plot(mean(temp), 'color', 'k', 'linewidth', 2)
                clusterEpilepsyRatioSA(i) = numElements / totalHoursSA;
                text(75, 20, strcat(num2str(clusterEpilepsyRatioSA(i)), {' events per hour'}));

            else
                patch([1:lengthLength, fliplr(1:lengthLength)], [patchVectorA, fliplr(patchVectorB)], 'r', 'FaceColor', 'r', 'EdgeColor', 'none', 'facealpha', .4) % , colorcolor, 'FaceColor', colorcolor, 'EdgeColor', 'none'
                plot(mean(temp), 'color', 'r', 'linewidth', 2)
                clusterEpilepsyRatioKA(i) = numElements / totalHoursKA;
                text(75, 20, strcat(num2str(clusterEpilepsyRatioKA(i)), {' events per hour'}));

            end
            xlim([0 lengthLength])
            ylim([-25 25])
    end
    suptitle(strcat({'Clustered Events from '}, treatmentGroup, {' Injected Animals'}))
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
    print(strcat('~/Desktop/clusterSummary_', treatmentGroup, '.pdf'), '-dpdf', '-painters')
    close all
end

clusterEpilepsyRatio_odds = clusterEpilepsyRatioKA ./ clusterEpilepsyRatioSA; % the odds ratio of being an event related to KA.
clusterEpilepsyRatio_prob = clusterEpilepsyRatio_odds ./ (clusterEpilepsyRatio_odds + 1) % 

% calculate and plot the hourly epilelspy index
h = figure('units', 'inch', 'pos', [10 10 5 5]) 

clear probsList
for i = 1:max(idx)
    probsList(idx == i) = clusterEpilepsyRatio_prob(i);
end

summaryTable = table(idx, probsList', eventLabel', eventTreatment');
summaryTable.Properties.VariableNames = {'idx', 'probsList', 'eventLabel', 'eventTreatment'};
summaryTable = grpstats(summaryTable, {'eventTreatment', 'eventLabel', 'idx'}, {'mean'}, 'DataVars',{'probsList'});
summaryTable.eventsPerHour = summaryTable.GroupCount;
for i = 1:size(summaryTable, 1)
    summaryTable.eventsPerHour(i) = summaryTable.GroupCount(i) / recordLengthsTable.sum_recordLengths(strcmp(summaryTable.eventLabel(i), recordLengthsTable.animalName));
end

summaryTable.mean_probsList = (summaryTable.mean_probsList - .5) .* 2;
summaryTable = summaryTable(summaryTable.mean_probsList > 0,:);
summaryTable.epilepsyEventsPerHour = summaryTable.eventsPerHour .* summaryTable.mean_probsList;

summaryTable = grpstats(summaryTable, {'eventTreatment', 'eventLabel'}, {'sum'}, 'DataVars',{'epilepsyEventsPerHour'});
summaryTable = summaryTable(~strcmp(summaryTable.eventTreatment, 'WT'), :);

[~, ppp,] = ttest2(summaryTable.sum_epilepsyEventsPerHour(strcmp(summaryTable.eventTreatment, 'KA')), summaryTable.sum_epilepsyEventsPerHour(strcmp(summaryTable.eventTreatment, 'SA')), 'Vartype','unequal')

boxplot(summaryTable.sum_epilepsyEventsPerHour, summaryTable.eventTreatment)
hold on
for i = 1:size(summaryTable.eventTreatment, 1)
    if strcmp(summaryTable.eventTreatment(i), 'KA')
        if strcmp(summaryTable.eventLabel(i), '1-5b_L1R0') || strcmp(summaryTable.eventLabel(i), '1-6b_L1R2')
            plot(normrnd(1, .05), summaryTable.sum_epilepsyEventsPerHour(i), 'm.', 'MarkerSize', 25)
        else
            plot(normrnd(1, .05), summaryTable.sum_epilepsyEventsPerHour(i), 'ro', 'MarkerSize', 7)
        end
    else
        plot(normrnd(2, .05), summaryTable.sum_epilepsyEventsPerHour(i), 'ko', 'MarkerSize', 7)
    end
end
ylabel('Hourly Epilepsy Index')
xlabel('Treatment')
text(1.5, 50, num2str(ppp))
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
print('~/Desktop/hourlyEpilepsyIndex.pdf', '-dpdf')

% plot examples of events from each cluster
for j = 1:2 
    if j == 1
        treatmentGroup = 'SA';
        treatmentColor = 'k';
    else
        treatmentGroup = 'KA';
        treatmentColor = 'r';
    end
    for i = 1:max(idx)
        h = figure('units', 'inch', 'pos', [10 10 15 15]) 
        temp = eventMatrix(intersect(find(idx == i), find(strcmp(eventTreatment, {treatmentGroup}))),:);
        nEvents = 100;
        if size(temp, 1) > nEvents
            randEvents = randsample(size(temp, 1), nEvents);
            temp = temp(randEvents,:);
        else
            nEvents = size(temp, 1);
        end
        for ii = 1:nEvents
            subplot(10, 10, ii)
                if j == 1
                    plot(temp(ii,:), 'k')
                else
                    plot(temp(ii,:), 'r')
                end
                %axis off
                xlim([0 size(eventMatrix, 2)])
                ylim([-20 20])

        end
        suptitle(strcat({'Examples of '}, treatmentGroup, {' events from Cluster: '}, num2str(i)))
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
        print(strcat('~/Desktop/examples_', treatmentGroup, '_cluster', num2str(i), '.pdf'), '-dpdf', '-painters')
        close all
    end
end
   
% Plot some of the features for each cluster -- 3d plot   
h = figure('units', 'inch', 'pos', [10 10 15 15]) 
for i = 1:max(idx)
    subplot(3, 3, i)
    
        temp2SA = eventAmplitude(intersect(find(idx == i), find(strcmp(eventTreatment, {'SA'}))));
        temp2KA = eventAmplitude(intersect(find(idx == i), find(strcmp(eventTreatment, {'KA'}))));
        temp3SA = eventDuration(intersect(find(idx == i), find(strcmp(eventTreatment, {'SA'}))));
        temp3KA = eventDuration(intersect(find(idx == i), find(strcmp(eventTreatment, {'KA'}))));
        scatter(temp2KA, temp3KA, '.r')
        hold on
        scatter(temp2SA, temp3SA, '.k')

        xlabel('Event Amplitdue')
        ylabel('Event Duration')
        xlim([5, 50])
        ylim([0, 1000])
        grid on
        title(strcat({'Cluster: '}, num2str(i)))
end
suptitle(strcat({'eventAmplitude by eventDuration'}))
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
print(strcat('~/Desktop/clusterSummary_duration_vs_amplitude.pdf'), '-dpdf', '-painters')

% boxplots that show the relative number of events per hour (figure 4B)
summaryTable = table(idx, probsList', eventLabel', eventTreatment');
summaryTable.Properties.VariableNames = {'idx', 'probsList', 'eventLabel', 'eventTreatment'};
summaryTable = grpstats(summaryTable, {'eventTreatment', 'eventLabel', 'idx', 'probsList'}, {'mean'}, 'DataVars',{'probsList'});
summaryTable.eventsPerHour = summaryTable.GroupCount;
for i = 1:size(summaryTable, 1)
    summaryTable.eventsPerHour(i) = summaryTable.GroupCount(i) / recordLengthsTable.sum_recordLengths(strcmp(summaryTable.eventLabel(i), recordLengthsTable.animalName));
end
summaryTable = summaryTable(~strcmp(summaryTable.eventTreatment, 'WT'), :);

% adds zeros to the summary Table.
uniqueAnimals = unique(summaryTable.eventLabel);
for i = 1:max(idx)
    for j = 1:length(uniqueAnimals)
       if sum(summaryTable.idx == i & strcmp(summaryTable.eventLabel, uniqueAnimals(j))) == 0
           tempEventTreatment = summaryTable.eventTreatment(strcmp(summaryTable.eventLabel, uniqueAnimals(j)));
           tempEventTreatment = tempEventTreatment(1);
           summaryTable = [summaryTable; {tempEventTreatment, uniqueAnimals(j), i, mean(summaryTable.probsList(summaryTable.idx == 1)), 0, 0, 0}]
       end
    end
end

    

h = figure('units', 'inch', 'pos', [10 10 15 15]) 
for i = 1:max(idx)
    subplot(3,3,i)
        summaryTableSubset = summaryTable(summaryTable.idx ==  i,:);
        boxplot(summaryTableSubset.eventsPerHour, summaryTableSubset.eventTreatment, 'color','k') 
        
        hold on
        [~, ppp,] = ttest2(summaryTableSubset.eventsPerHour(strcmp(summaryTableSubset.eventTreatment, 'KA')), summaryTableSubset.eventsPerHour(strcmp(summaryTableSubset.eventTreatment, 'SA')), 'Vartype','unequal')
        for zz = 1:size(summaryTableSubset.eventTreatment, 1)
            if strcmp(summaryTableSubset.eventTreatment(zz), 'KA')
                if strcmp(summaryTableSubset.eventLabel(zz), '1-5b_L1R0') || strcmp(summaryTableSubset.eventLabel(zz), '1-6b_L1R2')
                    plot(normrnd(1, .05), summaryTableSubset.eventsPerHour(zz), 'm.', 'MarkerSize', 25)
                else
                	plot(normrnd(1, .05), summaryTableSubset.eventsPerHour(zz), 'ro', 'MarkerSize', 7) 
                end
            else
                plot(normrnd(2, .05), summaryTableSubset.eventsPerHour(zz), 'ko', 'MarkerSize', 7)
            end
        end
        title(strcat({'Cluster Number: '}, num2str(i)))
        ylabel('Number of Events Per Hour')
        ylim([0 inf])
       
end

for i = 9
    subplot(3,3,i)
        summaryTableSubset = summaryTable(summaryTable.idx ==  i,:);
        boxplot(summaryTableSubset.eventsPerHour, summaryTableSubset.eventTreatment, 'color','k') 
        
        hold on
        [~, ppp,] = ttest2(summaryTableSubset.eventsPerHour(strcmp(summaryTableSubset.eventTreatment, 'KA')), summaryTableSubset.eventsPerHour(strcmp(summaryTableSubset.eventTreatment, 'SA')), 'Vartype','unequal')
        for zz = 1:size(summaryTableSubset.eventTreatment, 1)
            if strcmp(summaryTableSubset.eventTreatment(zz), 'KA')
                if strcmp(summaryTableSubset.eventLabel(zz), '1-5b_L1R0') || strcmp(summaryTableSubset.eventLabel(zz), '1-6b_L1R2')
                    plot(normrnd(1, .05), summaryTableSubset.eventsPerHour(zz), 'm.', 'MarkerSize', 25)
                else
                	plot(normrnd(1, .05), summaryTableSubset.eventsPerHour(zz), 'ro', 'MarkerSize', 7) 
                end
            else
                plot(normrnd(2, .05), summaryTableSubset.eventsPerHour(zz), 'ko', 'MarkerSize', 7)
            end
        end

        title(strcat({'Cluster Number: '}, num2str(i)))
        ylabel('Number of Events Per Hour')
        ylim([0 .06])
       
end
    
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
print(strcat('~/Desktop/anotherSummaryPlot_2.pdf'), '-dpdf', '-painters')
close all

% ----- Does sample rate make a difference? Nope.

% show some of the features and how they differ between the 256Hz and 512 Hz sample rates

% import namefile and bulid another dataFrame object to check the original sample rate
namefilespec = strcat(basePath, '/jonesLab_git/jonesLab_code/Projects/falconHawk/falconHawk_namefile_CURRENT.txt');
F = ParseNameFile(namefilespec); % makes the object 'F'

% build a dataFrame file including metadat
n = 1;
for i = 1:length(F.f)
    if str2num(F.itemparams{1, i}.analyzeEEG{1}) == 1
        % add edffilespecs for as many files as listed in the namefile
        for j = 1:length(F.itemparams{1, i}.EEGfile)
            dataFrame2(n).edffilespec = strcat(labDataDrive, '/jonesLab_data/falconHawk/eeg_data/', F.itemparams{1, i}.animalGroup{1} , filesep , F.itemparams{1, i}.EEGfile{1, j});
            [a, b, c] = fileparts(dataFrame2(n).edffilespec);
            dataFrame2(n).uniqueID = strcat(F.itemparams{1, i}.animalGroup{1}, '_', F.itemparams{1, i}.animalID{1});
            dataFrame2(n).fs = F.itemparams{1, i}.recordFs{1};
            dataFrame2(n).channel = F.itemparams{1, i}.EEGchan{1, j};
            dataFrame2(n).treatment = F.itemparams{1, i}.treatment{1};
            n = n + 1;
        end
    end
end

for i = 1:length(dataFrame2)
    which256(i) = strcmp(dataFrame2(i).fs, '256');
end
counter = 1;
for i = 1:length(dataFrame2)
    if which256(i) == 1
        nameWhich256{counter} = dataFrame2(i).uniqueID;
        counter = counter + 1;
    end
end
names256 = unique(nameWhich256);

summaryTable = table(idx, probsList', eventLabel', eventTreatment');
summaryTable.Properties.VariableNames = {'idx', 'probsList', 'eventLabel', 'eventTreatment'};
summaryTable = grpstats(summaryTable, {'eventTreatment', 'eventLabel', 'idx'}, {'mean'}, 'DataVars',{'probsList'});
summaryTable.eventsPerHour = summaryTable.GroupCount;
for i = 1:size(summaryTable, 1)
    summaryTable.eventsPerHour(i) = summaryTable.GroupCount(i) / recordLengthsTable.sum_recordLengths(strcmp(summaryTable.eventLabel(i), recordLengthsTable.animalName));
end
summaryTable = summaryTable(~strcmp(summaryTable.eventTreatment, 'WT'), :);


h = figure('units', 'inch', 'pos', [10 10 15 15]) 
for i = 1:max(idx)
    subplot(3,3,i)
        summaryTableSubset = summaryTable(summaryTable.idx ==  i,:);
        boxplot(summaryTableSubset.eventsPerHour, summaryTableSubset.eventTreatment, 'color','k') %
        
        hold on
        for zz = 1:size(summaryTableSubset.eventTreatment, 1)
            if strcmp(summaryTableSubset.eventTreatment(zz), 'KA')
                if sum(strcmp(summaryTableSubset.eventLabel(zz), names256)) > 0
                    plot(normrnd(1, .05), summaryTableSubset.eventsPerHour(zz), 'go', 'MarkerSize', 7)
                else
                	plot(normrnd(1, .05), summaryTableSubset.eventsPerHour(zz), 'ro', 'MarkerSize', 7) 
                end
            else
                if sum(strcmp(summaryTableSubset.eventLabel(zz), names256)) > 0
                    plot(normrnd(2, .05), summaryTableSubset.eventsPerHour(zz), 'go', 'MarkerSize', 7)
                else
                    plot(normrnd(2, .05), summaryTableSubset.eventsPerHour(zz), 'ko', 'MarkerSize', 7)
                end
            end
        end
        title(strcat({'Cluster Number: '}, num2str(i)))
        ylabel('Number of Events Per Hour')
        ylim([0 inf])
end
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  
print(strcat('~/Desktop/anotherSummaryPlot_2_256thing.pdf'), '-dpdf', '-painters')
close all


%% Find electrographic seizures

minSeizureLength = fs * 3; 
potentialSeizureList = find(eventDuration > minSeizureLength);

sum(strcmp(eventTreatment(potentialSeizureList), 'KA')) / 15
sum(strcmp(eventTreatment(potentialSeizureList), 'SA')) / 7

% show the eventMatrix snippets for each of the events identified above.
figure('units', 'inch', 'pos', [25 10 17 30]) 
counter = 1;
for i = potentialSeizureList
    subplot(length(potentialSeizureList), 1, counter)
        plot(eventMatrix(i,:))
    counter = counter + 1;
    axis off
end

% show the longer duratio clips.
figure('units', 'inch', 'pos', [25 10 17 30]) 
counter = 1;
for i = potentialSeizureList
    tempStart = eventStart(i);
    tempFilespec = eventFilespec(i);
    tempLabel = eventLabel(i);
    for j = 1:length(dataFrame)
        matchymatchy(j) = strcmp(tempFilespec, dataFrame(j).edffilespec);
        matchymatchy2(j) = strcmp(tempLabel, dataFrame(j).uniqueID);
    end
    theMatch = intersect(find(matchymatchy), find(matchymatchy2));
    subplot(length(potentialSeizureList), 1, counter)
        plot(dataFrame(theMatch).signal(tempStart:tempStart + (fs * 20)))
        axis off
        disp(dataFrame(theMatch).uniqueID)
    counter = counter + 1;
end

%%
% save model for future analyses
save('/Users/jesse/Google Drive/research/jonesLab/manuscripts/IIS detection manuscript/IISDetection_Share/Data from Jesse/IISFeatures_KAmodel.mat', 'PCA', 'eventDuration', 'clusterEpilepsyRatio_prob', 'gmBest')

    

%% save this environment?
%save(strcat(labDataDrive, '/jonesLab_data/falconHawk/eeg_data/currentWorkingFile.mat'), '-v7.3');
save('/Users/jesse/Google Drive/research/jonesLab/manuscripts/IIS detection manuscript/IISDetection_Share/Data from Jesse/falconHawk_EEGdata_512.mat', '-v7.3')


%% stuff not using anymore

%{
% Import the wavelet stuff from Rachel and make it into a matrix that we can append into the PCA data -- should eventually build this into the script here. It's only a few lines of code.
load('/Users/jesse/private/jonesLab_data/falconHawk/eeg_data/analysesFiles/wavedecout.mat')
counter = 1;
clear waveMat
waveMat = zeros(length(wavedecout), length(wavedecout(1).Cmatrix));
for i = 1:length(wavedecout)
    waveMat(counter,:) = wavedecout(i).Cmatrix;
    counter = counter + 1;
end

%save('~/Desktop/eventMatrix3_1_18.mat', 'eventMatrix', 'eventLabel', 'eventTreatment', 'eventDuration', 'eventAmplitude', 'waveMat') 
%}

