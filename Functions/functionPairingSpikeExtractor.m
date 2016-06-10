%This code will open matclust files, determine number of clusters, extract
%all spikes, and calculate RPVs, overall firing rate, and average
%waveforms. Overall, very similar to the functionSpikeWaveExtraction code,
%just different in terms of what it outputs.
function [masterStruct] = functionPairingSpikeExtractor(truncatedNames,...
    matclustFiles,rpvTime,clusterWindow,masterStruct);
%determine the number of clusters within each matclust file. Also extracts
%all spike times, rpvs, and waveforms.
numClusters = zeros(size(truncatedNames,2),1);
for i = 1:size(truncatedNames,2)
    %opens matclust file
    matclustFile = open(matclustFiles{i});
    %determines number of clusters within given file
    numClusters(i) = size(matclustFile.clustattrib.clustersOn,1); 
    %sets up empty arrays
    clusterSpikes = cell(numClusters(i),1);
    selectedSpikes = cell(numClusters(i),1);
    rpvViolationPercent = zeros(numClusters(i),1);
    %calculates inter-spike interval, saves this information
    for j = 1:numClusters(i)
        %this line pulls the actual indices of spikes for the cluster
        clusterSpikes{j} = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(j)}.index;
        %this subtracts all adjacent spikes
        diffSpikes = diff(matclustFile.clustdata.params(clusterSpikes{j},1));
        %this calculates all spikes within the refractory period violation
        %period
        rpvViolationPercent(j) = size(diffSpikes(diffSpikes<rpvTime),1)/size(clusterSpikes{j},1)*100;
        %this then windows out spikes, to only plot the small ISIs
        selectedSpikes{j} = diffSpikes(diffSpikes<clusterWindow(2)); %removes long pauses
        diffSpikes = [];
    end
    %prepares cluster indices. Finds actual points of index per cluster.
    %Then replaces with real times (in seconds)
    spikeTimes = cell(numClusters(i),1);
    for j = 1:numClusters(i)
        spikeTimes{j} = matclustFile.clustattrib.clusters{1,matclustFile.clustattrib.clustersOn(j)}.index;
        spikeTimes{j} = matclustFile.clustdata.params(spikeTimes{j},1);
    end
    %calculates overall firing rate
    totalTime = matclustFile.clustdata.datarange(2,1)-matclustFile.clustdata.datarange(1,1);
    overallFiring = zeros(numClusters(i),1);
    for j = 1:numClusters(i)
        overallFiring(j) = size(spikeTimes{j},1)/totalTime;
    end
    %pull out average waveform and standard error
    targetWaveName = strcat('waves',truncatedNames{i}(15:end),'.mat');
    waveLoader = open(targetWaveName);
    waveLoader =squeeze(waveLoader.waves);
    waveHolder = cell(numClusters(i),1);
    averageWaveHolder = zeros(size(waveLoader,1),numClusters(i),3);
    for j = 1:numClusters(i)
        waveHolder{j} = waveLoader(:,clusterSpikes{j});
        averageWaveHolder(:,j,2) = mean(waveHolder{j},2);
        averageWaveHolder(:,j,1) = mean(waveHolder{j},2)-std(waveHolder{j},0,2)/sqrt(size(waveHolder{j},2));
        averageWaveHolder(:,j,3) = mean(waveHolder{j},2)+std(waveHolder{j},0,2)/sqrt(size(waveHolder{j},2));
    end
    for j = 1:numClusters(i)
        waveHolder{j} = waveLoader(:,clusterSpikes{j});
        averageWaveHolder(:,j,2) = mean(waveHolder{j},2);
        averageWaveHolder(:,j,1) = mean(waveHolder{j},2)-std(waveHolder{j},0,2);
        averageWaveHolder(:,j,3) = mean(waveHolder{j},2)+std(waveHolder{j},0,2);
    end
    masterStruct.(truncatedNames{i}).Clusters = numClusters(i);
    masterStruct.(truncatedNames{i}).ISIData = selectedSpikes; 
    masterStruct.(truncatedNames{i}).RPVs = rpvViolationPercent;
    masterStruct.(truncatedNames{i}).SpikeTimes = spikeTimes;
    masterStruct.(truncatedNames{i}).OverallFiringRate = overallFiring;
    masterStruct.(truncatedNames{i}).AverageWaveForms = averageWaveHolder;
    masterStruct.(truncatedNames{i}).AllWaveForms = waveHolder;
end

end