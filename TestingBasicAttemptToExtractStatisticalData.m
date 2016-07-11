holder = matclustStruct;

targetCellRasters = holder.matclust_param_nt13.Rasters{1,1};

rasterWindow = [-0.5,0.5];
histBin = 0.01; %bin size in seconds
histBinNum = (rasterWindow(2)-rasterWindow(1))/histBin;
histBinVector = [rasterWindow(1)+histBin/2:histBin:rasterWindow(2)-histBin/2]; %this is vector with midpoints of all histogram bins

freqSize = size(holder.UniqueFreqs,1);
dbSize = size(holder.UniqueDBs,1);
repSize = holder.ToneReps;

responseCellArray = cell(6,13);
averageResp = zeros(6,13);
averageSTD = zeros(6,13);

indivTrialRasters = cell(6,13);
trialHolder = cell(repSize,1);

for i = 1:freqSize
    for j = 1:dbSize
        %the code in this part is to extract responses within the timing
        %window, and then calculate all responses within the tone
        %presentation period, so I can calculate average spikes and
        %standard deviation.
        trialNumHolder = intersect(find(holder.Frequencies == holder.UniqueFreqs(i)),...
            find(holder.dBs == holder.UniqueDBs(j)));
        responseHolder = zeros(size(trialNumHolder,1),2);
        responseHolder(:,1) = trialNumHolder;
        for k = 1:repSize
            indivResponse = targetCellRasters(targetCellRasters(:,3)==holder.UniqueFreqs(i) & ...
            targetCellRasters(:,4)==holder.UniqueDBs(j) &...
            targetCellRasters(:,1)==responseHolder(k,1),2);
        %This pulls out total number of spikes in the response period
            responseHolder(k,2) = size(find(indivResponse > 0 & indivResponse < holder.ToneDur),1);
            %This pulls out the actual spike times for calculating errors. 
            trialHolder{k} = indivResponse;
        end
        averageResp(j,i) = mean(responseHolder(:,2));
        averageSTD(j,i) = std(responseHolder(:,2));
        responseCellArray{j,i} = responseHolder;
        indivTrialRasters{j,i} = trialHolder;        
    end
end


%below is code for comparing different responses (each individual response)
%between two conditions. 
statsHolder = zeros(6,13);
for i = 1:6
    for j = 1:13
        statsHolder(i,j) = signrank(adjustedResponses{i,j}(:,2),origResponses{i,j}(:,2))
    end
end