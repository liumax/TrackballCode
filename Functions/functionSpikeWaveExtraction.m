%This code will be for extraction of spikes, calculation of RPVs, and
%calculation of overall firing rate. This code will also extract the
%average waveform of the cluster.

function [matclustStruct, clusterSizer] = functionSpikeWaveExtraction(rpvTime,...
    i,matclustFiles,matclustStruct,truncatedNames,clusterWindow);

%Opens the matclust file.
matclustFile = open(matclustFiles{i});
%Now I need to actually pull components of the clusters:
clusterSizer= size(matclustFile.clustattrib.clustersOn,1); %pulls number of clusters

    
%preps series of arrays for extracting spikes
clusterSpikes = cell(clusterSizer,1);
selectedSpikes = cell(clusterSizer,1);
rpvViolationPercent = zeros(clusterSizer,1);
rpvNumber = zeros(clusterSizer,1);
totalNumber = zeros(clusterSizer,1);

%calculates inter-spike interval, saves this information
for j = 1:clusterSizer
    %this line pulls the actual indices of spikes for the cluster
    clusterSpikes{j} = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(j)}.index;
    %this subtracts all adjacent spikes
    diffSpikes = diff(matclustFile.clustdata.params(clusterSpikes{j},1));
    %this calculates all spikes within the refractory period violation
    %period
    rpvNumber(j) = size(diffSpikes(diffSpikes<rpvTime),1);
    totalNumber(j) = size(clusterSpikes{j},1);
    rpvViolationPercent(j) = size(diffSpikes(diffSpikes<rpvTime),1)/size(clusterSpikes{j},1)*100;
    %this then windows out spikes, to only plot the small ISIs
    selectedSpikes{j} = diffSpikes(diffSpikes<clusterWindow(2)); %removes long pauses
    diffSpikes = [];
end

%prepares cluster indices. Finds actual points of index per cluster.
%Then replaces with real times (in Trodes time)
spikeTimes = cell(clusterSizer,1);
for j = 1:clusterSizer
    spikeTimes{j} = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(j)}.index;
    spikeTimes{j} = matclustFile.clustdata.params(spikeTimes{j},1);
end

%calculates overall firing rate
totalTime = matclustFile.clustdata.datarange(2,1)-matclustFile.clustdata.datarange(1,1);
overallFiring = zeros(clusterSizer,1);
for j = 1:clusterSizer
    overallFiring(j) = size(spikeTimes{j},1)/totalTime;
end

%pull out average waveform and standard error
targetWaveName = strcat('waves_',truncatedNames{i},'.mat');
waveLoader = open(targetWaveName);
waveLoader =squeeze(waveLoader.waves);
waveHolder = cell(clusterSizer,1);
averageWaveHolder = zeros(size(waveLoader,1),clusterSizer,3);
for j = 1:clusterSizer
    waveHolder{j} = waveLoader(:,clusterSpikes{j});
    averageWaveHolder(:,j,2) = mean(waveHolder{j},2);
    averageWaveHolder(:,j,1) = mean(waveHolder{j},2)-std(waveHolder{j},0,2)/sqrt(size(waveHolder{j},2));
    averageWaveHolder(:,j,3) = mean(waveHolder{j},2)+std(waveHolder{j},0,2)/sqrt(size(waveHolder{j},2));
end
for j = 1:clusterSizer
    waveHolder{j} = waveLoader(:,clusterSpikes{j});
    averageWaveHolder(:,j,2) = mean(waveHolder{j},2);
    averageWaveHolder(:,j,1) = mean(waveHolder{j},2)-std(waveHolder{j},0,2);
    averageWaveHolder(:,j,3) = mean(waveHolder{j},2)+std(waveHolder{j},0,2);
end

matclustStruct.(truncatedNames{i}).Clusters =  clusterSizer; %stores this number into structured array
matclustStruct.(truncatedNames{i}).ISIData = selectedSpikes;
matclustStruct.(truncatedNames{i}).RPVs = rpvViolationPercent;
matclustStruct.(truncatedNames{i}).RPVNumber = rpvNumber;
matclustStruct.(truncatedNames{i}).TotalSpikeNumber = totalNumber;
matclustStruct.(truncatedNames{i}).SpikeTimes = spikeTimes;
matclustStruct.(truncatedNames{i}).OverallFiringRate = overallFiring;
matclustStruct.(truncatedNames{i}).AverageWaveForms = averageWaveHolder;
matclustStruct.(truncatedNames{i}).AllWaveForms = waveHolder;
matclustStruct.(truncatedNames{i}).TotalTimeRecording = totalTime;
end

















