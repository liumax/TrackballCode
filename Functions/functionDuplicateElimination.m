



function [s] = functionDuplicateElimination(s);

%pull parameters
downSampFactor = s.Parameters.DownSampFactor;
corrSlide = s.Parameters.corrSlide;
threshComp = s.Parameters.ThresholdComparison;

%pull time information to construct 
firstPoint = s.TimeFilterRange(1);
lastPoint = s.TimeFilterRange(2);

%compute lag window for xcorr
lagWindow = corrSlide*s.Parameters.trodesFS/downSampFactor;
lags = [-lagWindow:1:lagWindow];

%pull number of units
numUnits = length(s.DesignationArray);
names = s.DesignationName;

%construct spike trains, store in temporary structure
tmp = struct;
for tmpCount = 1:numUnits
    tmp.(names{tmpCount}).SpikeTrain = zeros(round((lastPoint-firstPoint)*30000/downSampFactor),1);
    tmp.(names{tmpCount}).SpikeTrain(round((s.(names{tmpCount}).SpikeTimes-firstPoint)*30000/downSampFactor)) = 1;
    %also create a temporally downsampled list of spike times
    tmp.(names{tmpCount}).SpikeTimes = round((s.(names{tmpCount}).SpikeTimes-firstPoint)*30000/downSampFactor);
end

%figure out what kind of probe I'm using!
numTrodes = s.NumberTrodes;
if numTrodes == 4;
    %this means i'm on the 16 channel single shank
    shanks = 1;
elseif numTrodes == 8;
    %this means i'm on the 32 channel double shank
    shanks = 2;
elseif numTrodes == 16;
    %this means i'm on the 64 channel double shank
    shanks = 2;
end

%this generates an array for the shanks.
shankDesig = zeros(numTrodes/shanks,shanks);
holderShank = 1;
for shnkInd = 1:shanks
    shankDesig(:,shnkInd) = [holderShank:holderShank+numTrodes/shanks-1];
    holderShank = holderShank + numTrodes/shanks;
end

%set subplot settings!
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);

%now go through systematically and look at cross correlations. What we are
%doing is advancing forward down the shank, never backwards (with the
%assumption that any correlations were already taken care of). 

%pull figure. 
hFig = figure;
set(hFig, 'Position', [10 80 600 600])

for shnkInd = 1:shanks
    disp(strcat('Analyzing Shank:',num2str(shnkInd)))
    %find the correct names for the targeted shank
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
                disp(strcat('Analyzing:',names{uniqueUnits(unitInd)}))
                %first step is gross go over, will basically use the crude
                %old code I was using to see if there is any chance of
                %overlap between the two units. 
                compHolder = zeros(length(indNeighbor),3);
                for nghbrInd = 1:length(indNeighbor)
                    compHolder(nghbrInd,1) = length(intersect(tmp.(names{uniqueUnits(unitInd)}).SpikeTimes,tmp.(names{indNeighbor(nghbrInd)}).SpikeTimes));
                    compHolder(nghbrInd,2) = compHolder(nghbrInd,1)/length(tmp.(names{uniqueUnits(unitInd)}).SpikeTimes);
                    compHolder(nghbrInd,3) = compHolder(nghbrInd,1)/length(tmp.(names{indNeighbor(nghbrInd)}).SpikeTimes);
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
                    disp(strcat(num2str(length(compFinder)),' Potential Duplicates Detected, Performing xCorr'))
                    while ~isempty(compFinder)
                        %computer xcorr
                        [r,~] = xcorr(tmp.(names{uniqueUnits(unitInd)}).SpikeTrain,...
                            tmp.(names{indNeighbor(compFinder(1))}).SpikeTrain,lagWindow);
                        %plot to display
                        clf;
                        subplot(3,1,1)
                        plot(lags,r)
                        xlim([-lagWindow,lagWindow])
                        ylim([0 max(r)])
                        title(strcat('Correlogram of:',names{uniqueUnits(unitInd)},'&',names{indNeighbor(compFinder(1))}))
                        subplot(3,1,2)
                        hold on
                        trueMax = 0;
                        trueMin = 0;
                        numWaves = size(s.(names{uniqueUnits(unitInd)}).AverageWaveForms,2);
                        for m =1:numWaves
                            plot(s.(names{uniqueUnits(unitInd)}).AverageWaveForms(:,m)-(100*(m-1)),'LineWidth',2)
                            testMax = max(s.(names{uniqueUnits(unitInd)}).AverageWaveForms(:,m)-(100*(m-1)));
                            if testMax > trueMax
                                trueMax = testMax;
                            end
                            testMin = min(s.(names{uniqueUnits(unitInd)}).AverageWaveForms(:,m)-(100*(m-1)));
                            if testMin < trueMin
                                trueMin = testMin;
                            end
                        end
                        xlim([0 length(s.(names{uniqueUnits(unitInd)}).AverageWaveForms(:,1))])
                        ylim([trueMin trueMax])
                        title(strcat('Waves of:',names{uniqueUnits(unitInd)}))
                        subplot(3,1,3)
                        hold on
                        trueMax = 0;
                        trueMin = 0;
                        for m =1:numWaves
                            plot(s.(names{indNeighbor(compFinder(1))}).AverageWaveForms(:,m)-(100*(m-1)),'LineWidth',2)
                            testMax = max(s.(names{indNeighbor(compFinder(1))}).AverageWaveForms(:,m)-(100*(m-1)));
                            if testMax > trueMax
                                trueMax = testMax;
                            end
                            testMin = min(s.(names{indNeighbor(compFinder(1))}).AverageWaveForms(:,m)-(100*(m-1)));
                            if testMin < trueMin
                                trueMin = testMin;
                            end
                        end
                        xlim([0 length(s.(names{indNeighbor(compFinder(1))}).AverageWaveForms(:,1))])
                        ylim([trueMin trueMax])
                        title(strcat('Waves of:',names{indNeighbor(compFinder(1))}))
                        %% Ask for input
                        promptCounter = 1; %This is used to run the while loop.
                        whileCounter = 0; %this is the counter that gets updated to exit the loop

                        while whileCounter ~= promptCounter
                            try
                                prompt = 'Decide on Units: k: keep both f: keep first s: keep second';
                                str = input(prompt,'s');
                                if str~='k' & str~='f' & str~='s'
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
                            targetUnits(1)=[];
                            unitInd = unitInd + 1;
                        elseif strfind(str,'f')
                            disp('Keeping First, Deleting Second')
                            %create new variable in s to record unit as
                            %being deleted. Save the unit that was the
                            %reason for deletion. This is to say the field
                            %will be the deleted unit. The name in the
                            %field will be the unit that triggered the
                            %deletion. 
                            s.DeletedUnits.(names{indNeighbor(compFinder(1))}) = (names{uniqueUnits(unitInd)});
                            %delete data from second unit from s
                            s = rmfield(s,(names{indNeighbor(compFinder(1))}));
                            %remove from s names and designation array
                            %find the target name
                            nameInd = find(~cellfun(@isempty,strfind(s.DesignationName,names{indNeighbor(compFinder(1))})));
                            s.DesignationName(nameInd) = [];
                            names(nameInd) = [];
                            s.DesignationArray(nameInd,:) = [];
                            %increment compFinder
                            compFinder(1) = [];
                            targetUnits(1)=[];
                            unitInd = unitInd + 1;
                        elseif strfind(str,'s')
                            disp('Keeping Second, Deleting First')
                            %create new variable in s to record unit as
                            %being deleted. Save the unit that was the
                            %reason for deletion
                            s.DeletedUnits.(names{uniqueUnits(unitInd)}) = (names{indNeighbor(compFinder(1))});
                            %delete data from second unit from s
                            s = rmfield(s,(names{uniqueUnits(unitInd)}));
                            %remove from s names and designation array
                            nameInd = find(~cellfun(@isempty,strfind(s.DesignationName,names{uniqueUnits(unitInd)})));
                            s.DesignationName(nameInd) = [];
                            names(nameInd) = [];
                            s.DesignationArray(nameInd,:) = [];
                            %delete compFinder, since the target of
                            %comparison has been removed
                            compFinder = [];
                            targetUnits(1)=[];
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

                
                












