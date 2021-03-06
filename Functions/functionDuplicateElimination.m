



function [s] = functionDuplicateElimination(s,downSampFactor,corrSlide,threshComp,...
    samplingRate,rpvTime,clusterWindow,shankDesig,shankMap,shanks);

%pull parameters
% downSampFactor = s.Parameters.DownSampFactor;
% corrSlide = s.Parameters.corrSlide;
% threshComp = s.Parameters.ThresholdComparison;

%pull time information to construct 
firstPoint = s.TimeFilterRange(1);
lastPoint = s.TimeFilterRange(2);

%compute lag window for xcorr
lagWindow = corrSlide;
crossWindow = [-lagWindow,lagWindow];
lags = [-lagWindow:0.001:lagWindow];

%pull number of units
numUnits = length(s.DesignationArray);

%construct spike trains, store in temporary structure
tmp = struct;
for tmpCount = 1:numUnits
    %also create a temporally downsampled list of spike times
    tmp.(s.DesignationName{tmpCount}).SpikeTimes = round((s.(s.DesignationName{tmpCount}).SpikeTimes-firstPoint)*30000/downSampFactor);
end

%set subplot settings!
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.04], [0.05 0.08], [0.06 0.03]);

%now go through systematically and look at cross correlations. What we are
%doing is advancing forward down the shank, never backwards (with the
%assumption that any correlations were already taken care of). 

%pull figure. 
hFig = figure;
set(hFig, 'Position', [10 80 600 600])

for shnkInd = 1:shanks
    disp(strcat('Analyzing Shank:',num2str(shnkInd)))
    %find the correct s.DesignationName for the targeted shank
    indShank = find(ismember(s.DesignationArray(:,1),shankDesig(:,shnkInd)));
    %find the unique trodes represented
    uniqueTrodes = unique(s.DesignationArray(indShank,1));
    for trdInd = 1:length(uniqueTrodes)
        disp(strcat('Analyzing Trode:',num2str(trdInd)))
        %pulls unit identifiers
        uniqueUnits = find(s.DesignationArray(:,1) == uniqueTrodes(trdInd));
        %find if there are nearby tetrodes. Only go one above
        neighbor = uniqueTrodes(trdInd)+1;
        %check if this trode is represented
        neighborCheck = double(ismember(neighbor,uniqueTrodes));
        if neighborCheck == 1
            %now we go through each unit of the current array to compare it
            %with neighbors
            targetUnits = uniqueUnits;
            unitInd = 1;
            while ~isempty(targetUnits)
                %find the indices for these units
                indNeighbor = find(ismember(s.DesignationArray(:,1),neighbor));
                uniqueUnits = find(s.DesignationArray(:,1) == uniqueTrodes(trdInd));
                disp(strcat('Analyzing:',s.DesignationName{uniqueUnits(unitInd)}))
                %first step is gross go over, will basically use the crude
                %old code I was using to see if there is any chance of
                %overlap between the two units. 
                compHolder = zeros(length(indNeighbor),3);
                for nghbrInd = 1:length(indNeighbor)
                    compHolder(nghbrInd,1) = length(intersect(tmp.(s.DesignationName{uniqueUnits(unitInd)}).SpikeTimes,tmp.(s.DesignationName{indNeighbor(nghbrInd)}).SpikeTimes));
                    compHolder(nghbrInd,2) = compHolder(nghbrInd,1)/length(tmp.(s.DesignationName{uniqueUnits(unitInd)}).SpikeTimes);
                    compHolder(nghbrInd,3) = compHolder(nghbrInd,1)/length(tmp.(s.DesignationName{indNeighbor(nghbrInd)}).SpikeTimes);
                end
                %find maximum percentage of spikes shared.
                compTester = max(compHolder(:,[2:3]),[],2);
                %find any that exceed threshold
                compFinder = find(compTester > threshComp);
                if isempty(compFinder)
                    disp('No Duplicates Detected')
                    targetUnits(1)=[];
                    unitInd = unitInd + 1;
                elseif ~isempty(compFinder)
                    disp(strcat(num2str(length(compFinder)),' Potential Duplicates Detected, Performing OwnCorr'))
                    while ~isempty(compFinder)
                        %compute xcorr
                        spikeTimes1 = tmp.(s.DesignationName{uniqueUnits(unitInd)}).SpikeTimes;
                        spikeTimes2 = tmp.(s.DesignationName{indNeighbor(compFinder(1))}).SpikeTimes;
                        histStore = ones(100000,1);
                        histCounter = 1;
                        for spikeInd = 1:length(spikeTimes1)
                            %subtract spike time from first train out of all of second train
                            subSpikes = spikeTimes2 - spikeTimes1(spikeInd);
                            %remove things outside the window of interest
                            subSpikes(subSpikes>crossWindow(2) | subSpikes<crossWindow(1)) = [];
                            histStore(histCounter:(histCounter + length(subSpikes)-1)) = subSpikes;
                            histCounter = histCounter + length(subSpikes);
                        end
                        histStore(histCounter:end) = [];
%                         [r,~] = xcorr(tmp.(s.DesignationName{uniqueUnits(unitInd)}).SpikeTrain,...
%                             tmp.(s.DesignationName{indNeighbor(compFinder(1))}).SpikeTrain,lagWindow);
                        %% plot to display
                        clf;
                        subplot(3,1,1)
                        hist(histStore,lags)
                        xlim([-lagWindow,lagWindow])
                        title({strcat('Correlogram of:',s.DesignationName{uniqueUnits(unitInd)},'&',s.DesignationName{indNeighbor(compFinder(1))});...
                            strcat(s.DesignationName{uniqueUnits(unitInd)},':',num2str(length(tmp.(s.DesignationName{uniqueUnits(unitInd)}).SpikeTimes)),',',...
                            s.DesignationName{indNeighbor(compFinder(1))},':',num2str(length(tmp.(s.DesignationName{indNeighbor(compFinder(1))}).SpikeTimes)))})
                        subplot(3,2,3)
                        hold on
                        trueMax = 0;
                        trueMin = 0;
                        numWaves = size(s.(s.DesignationName{uniqueUnits(unitInd)}).AverageWaveForms,2);
                        for m =1:numWaves
                            plot(s.(s.DesignationName{uniqueUnits(unitInd)}).AverageWaveForms(:,m)-(100*(m-1)),'LineWidth',2)
                            testMax = max(s.(s.DesignationName{uniqueUnits(unitInd)}).AverageWaveForms(:,m)-(100*(m-1)));
                            if testMax > trueMax
                                trueMax = testMax;
                            end
                            testMin = min(s.(s.DesignationName{uniqueUnits(unitInd)}).AverageWaveForms(:,m)-(100*(m-1)));
                            if testMin < trueMin
                                trueMin = testMin;
                            end
                        end
                        xlim([0 length(s.(s.DesignationName{uniqueUnits(unitInd)}).AverageWaveForms(:,1))])
                        ylim([trueMin trueMax])
                        title(strcat('Waves of:',s.DesignationName{uniqueUnits(unitInd)}))
                        
                        subplot(3,2,4)
                        hist(s.(s.DesignationName{uniqueUnits(unitInd)}).ISIGraph,1000)
                        histMax = max(hist(s.(s.DesignationName{uniqueUnits(unitInd)}).ISIGraph,1000));
                        line([rpvTime rpvTime],[0 histMax],'LineWidth',1,'Color','red')
                        xlim(clusterWindow)
                        title('First Unit AutoCorr')
                        
                        subplot(3,2,5)
                        hold on
                        trueMax = 0;
                        trueMin = 0;
                        for m =1:numWaves
                            plot(s.(s.DesignationName{indNeighbor(compFinder(1))}).AverageWaveForms(:,m)-(100*(m-1)),'LineWidth',2)
                            testMax = max(s.(s.DesignationName{indNeighbor(compFinder(1))}).AverageWaveForms(:,m)-(100*(m-1)));
                            if testMax > trueMax
                                trueMax = testMax;
                            end
                            testMin = min(s.(s.DesignationName{indNeighbor(compFinder(1))}).AverageWaveForms(:,m)-(100*(m-1)));
                            if testMin < trueMin
                                trueMin = testMin;
                            end
                        end
                        xlim([0 length(s.(s.DesignationName{indNeighbor(compFinder(1))}).AverageWaveForms(:,1))])
                        ylim([trueMin trueMax])
                        title(strcat('Waves of:',s.DesignationName{indNeighbor(compFinder(1))}))
                        
                        subplot(3,2,6)
                        hist(s.(s.DesignationName{indNeighbor(compFinder(1))}).ISIGraph,1000)
                        histMax = max(hist(s.(s.DesignationName{indNeighbor(compFinder(1))}).ISIGraph,1000));
                        line([rpvTime rpvTime],[0 histMax],'LineWidth',1,'Color','red')
                        xlim(clusterWindow)
                        title('First Unit AutoCorr')
                        %% Ask for input
                        promptCounter = 1; %This is used to run the while loop.
                        whileCounter = 0; %this is the counter that gets updated to exit the loop

                        while whileCounter ~= promptCounter
                            try
                                prompt = 'Decide on Units: k: keep both f: keep first s: keep second n: discard both';
                                str = input(prompt,'s');
                                if str~='k' & str~='f' & str~='s' & str~='n'
                                    error
                                else
                                    whileCounter = 1;
                                end
                            catch
                            end
                        end
                        if strfind(str,'k')
                            disp('Keep Both Units')
                            compFinder(1) = [];
                            %if statement is to protect against cases with
                            %multiple comparisons.
                            if isempty(compFinder)
                                targetUnits(1)=[];
                                unitInd = unitInd + 1;
                            end
                        elseif strfind(str,'f')
                            disp('Keeping First, Deleting Second')
                            %create new variable in s to record unit as
                            %being deleted. Save the unit that was the
                            %reason for deletion. This is to say the field
                            %will be the deleted unit. The name in the
                            %field will be the unit that triggered the
                            %deletion. 
                            s.DeletedUnits.(s.DesignationName{indNeighbor(compFinder(1))}) = (s.DesignationName{uniqueUnits(unitInd)});
                            %delete data from second unit from s
                            s = rmfield(s,(s.DesignationName{indNeighbor(compFinder(1))}));
                            %remove from s s.DesignationName and designation array
                            %find the target name
                            nameInd = find(~cellfun(@isempty,strfind(s.DesignationName,s.DesignationName{indNeighbor(compFinder(1))})));
                            s.DesignationName(nameInd) = [];
%                             s.DesignationName(nameInd) = [];
                            s.DesignationArray(nameInd,:) = [];
                            %increment compFinder
                            compFinder(1) = [];
                            %adust indNeighbor
                            indNeighbor = indNeighbor - 1;
                            %if statement is to protect against cases with
                            %multiple comparisons.
                            if isempty(compFinder)
                                targetUnits(1)=[];
                                unitInd = unitInd + 1;
                            end
                        elseif strfind(str,'s')
                            disp('Keeping Second, Deleting First')
                            %create new variable in s to record unit as
                            %being deleted. Save the unit that was the
                            %reason for deletion
                            s.DeletedUnits.(s.DesignationName{uniqueUnits(unitInd)}) = (s.DesignationName{indNeighbor(compFinder(1))});
                            %delete data from second unit from s
                            s = rmfield(s,(s.DesignationName{uniqueUnits(unitInd)}));
                            %remove from s s.DesignationName and designation array
                            nameInd = find(~cellfun(@isempty,strfind(s.DesignationName,s.DesignationName{uniqueUnits(unitInd)})));
                            s.DesignationName(nameInd) = [];
%                             s.DesignationName(nameInd) = [];
                            s.DesignationArray(nameInd,:) = [];
                            %delete compFinder, since the target of
                            %comparison has been removed
                            compFinder = [];
                            targetUnits(1)=[];
                            %dont need to advance unitInd since I'm
                            %deleting 
                        elseif strfind(str,'n')
                            disp('Deleting Both Units')
                            %create new variable in s to record units as
                            %being deleted. 
                            s.DeletedUnits.(s.DesignationName{uniqueUnits(unitInd)}) = (s.DesignationName{indNeighbor(compFinder(1))});
                            s.DeletedUnits.(s.DesignationName{indNeighbor(compFinder(1))}) = (s.DesignationName{uniqueUnits(unitInd)});
                            %deleted data from targeted units
                            s = rmfield(s,(s.DesignationName{indNeighbor(compFinder(1))}));
                            s = rmfield(s,(s.DesignationName{uniqueUnits(unitInd)}));
                            
                            %remove from s s.DesignationName and designation array
                            %find the target name
                            nameInd = find(~cellfun(@isempty,strfind(s.DesignationName,s.DesignationName{indNeighbor(compFinder(1))})));
                            s.DesignationName(nameInd) = [];
%                             s.DesignationName(nameInd) = [];
                            s.DesignationArray(nameInd,:) = [];
                            
                            %remove from s s.DesignationName and designation array
                            nameInd = find(~cellfun(@isempty,strfind(s.DesignationName,s.DesignationName{uniqueUnits(unitInd)})));
                            s.DesignationName(nameInd) = [];
%                             s.DesignationName(nameInd) = [];
                            s.DesignationArray(nameInd,:) = [];
                            %delete compFinder, since the target of
                            %comparison has been removed
                            compFinder = [];
                            targetUnits(1)=[];
                            %dont need to advance unitInd since I'm
                            %deleting 
                            
                            
                        end
                    end
                end
            end
        else
            disp('No Recorded Neighbors')
        end
    end
end


end

                
                












