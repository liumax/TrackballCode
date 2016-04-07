holder = matclustStruct;

targetCellRasters = holder.matclust_param_nt15.Rasters{1,1};

freqSize = size(holder.UniqueFreqs,1);
dbSize = size(holder.UniqueDBs,1);
repSize = holder.ToneReps;

responseCellArray = cell(6,13);
averageResp = zeros(6,13);
averageSTD = zeros(6,13);

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
            indivResponse = size(find(targetCellRasters(:,2)>0 & ...
            targetCellRasters(:,2)<holder.ToneDur & ...
            targetCellRasters(:,3)==holder.UniqueFreqs(i) & ...
            targetCellRasters(:,4)==holder.UniqueDBs(j) &...
            targetCellRasters(:,1)==responseHolder(k,1)),1);
            responseHolder(k,2) = indivResponse;
        end
        averageResp(j,i) = mean(responseHolder(:,2))/repSize;
        averageSTD(j,i) = std(responseHolder(:,2));
        responseCellArray{j,i} = responseHolder;
        
        %now I want to extract individual traces for frequency/amplitude. I
        %will do this by extracting all points from rasters and making them
        %into histograms. 
        extractHolder = targetCellRasters(targetCellRasters(:,3) == holder.UniqueFreqs(i) &...
            targetCellRasters(:,4) == holder.UniqueDBs(j),2);
        
        hist(extractHolder
        
    end
end
