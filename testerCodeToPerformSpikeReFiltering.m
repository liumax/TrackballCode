%This is code that is meant to re-pick out waveforms based on spike
%matching so that I can re-filter existing data. This code is meant to be
%used with the tuning dataset
%(FullTuningAnalysis170420RemovedCellsWithRepeats)

%This code must start in the mother directory with the folders of the data
%I want to analyze

load('FullTuningAnalysis170420RemovedCellsWithRepeats.mat')
fieldNames = fields(fullMaster);

spikeLim = [-0.0009 0.0009];

for i = 1:length(fieldNames)
    recordData = fullMaster.(fieldNames{i});
    recordName = fieldNames{i}(2:end) %this is the name for the file
    folderName = recordName(1:end-2); %this will be the name of the target folder
    mFolderName = strcat(recordName,'.matclust'); %this is the name for the targeted matclust folder name
    numUnits = length(recordData.DesignationName);
    desigNames = recordData.DesignationName;
    for j = 1:length(desigNames)
        disp(desigNames{j})
        unitData = recordData.(desigNames{j}).SpikeTimes;
        clustFind = strfind(desigNames{j},'cluster');
        paramName = strcat('param_',desigNames{j}(1:clustFind-1),'.mat');
        waveName = strcat('waves_',desigNames{j}(1:clustFind-1),'.mat');
        %move to the targeted folder
        cd (folderName)
        cd (mFolderName)
        %open target params file
        load((paramName))
        paramSpikes = filedata.params(:,1); %this should be the time column
        clear 'filedata' %close filedata
        %first, lets use intersect to try and pull out exactly matching
        %times
        waveWarn = ones(length(unitData),1); %this will be to warn about waves that arent perfect matches
        [C ia ib] = intersect(paramSpikes,unitData);
        %lets check to see if this covers all spikes
        if length(C) == length(unitData);
            %this is the case in which all the spikes have been covered. In
            %this case, we simply want to extract ia, which is the index of
            %total spike listing. This index will give us our desired
            %waveforms. 
            waveWarn(1:end) = 0;
            
            load((waveName))
            fullMaster.(fieldNames{i}).(desigNames{j}).OldWaves = fullMaster.(fieldNames{i}).(desigNames{j}).AllWaves;
            fullMaster.(fieldNames{i}).(desigNames{j}).AllWaves = [];
            fullMaster.(fieldNames{i}).(desigNames{j}).AllWaves = waves(:,:,ia);
            fullMaster.(fieldNames{i}).(desigNames{j}).AverageWaveForms = [];
            fullMaster.(fieldNames{i}).(desigNames{j}).AverageWaveForms = mean(waves(:,:,ia),3);
            fullMaster.(fieldNames{i}).(desigNames{j}).WaveWarning = waveWarn;
            clear 'waves'
            
        else %this is the case in which not all spikes are perfectly matched
            disp(strcat('Imperfect Match:',num2str(length(C)/length(unitData)*100),'%,',num2str(length(C)-length(unitData)),' Points Remaining'))
            %make list of times of remaining waves
            remainWaves = unitData;
            remainWaves(ib) = [];
            %adjust waveWarn
            waveWarn(ib) = 0;
            findWarn = find(waveWarn == 1);
            for k = 1:length(remainWaves)
                newLim = spikeLim; %use this as a placeholder so I can adjust the limits
                whileTrigger = 1
                while whileTrigger = 1;
                    waveDiff = paramSpikes - remainWaves(k);
                    diffFind = find(waveDiff > newLim(1) & waveDiff < newLim(2));
                    if isempty(diffFind) %if no spikes found, want to expand window
                        newLim(1) = newLim(1) - 0.0001;
                        newLim(2) = newLim(2) + 0.0001;
                    else %in the case where you do find spikes in the window
                        if length(diffFind == 1) %only one spike. Accept this spike
                            waveWarn(findWarn(k)) = waveDiff(diffFind); %update waveWarn
                            ia(end+1) = diffFind;
                            newLim = spikeLim;  %reset newLim
                            whileTrigger = 0; %exit while loop
                        elseif length(diffFind > 1) %more than one spike in the window! take the closest spike
                            absDiff = abs(waveDiff(diffFind));
                            [minVal minInd] = min(absDiff);
                            waveWarn(findWarn(k)) = waveDiff(diffFind(minInd)); %update waveWarn
                            ia(end+1) = diffFind(minInd);
                            newLim = spikeLim; %reset newLim
                            whileTrigger = 0; %exit while loop
                        end
                            
                    end         
                end
            end
            %now we've completed the finding of the waves! yay. 
            disp('Remaining Waves Found')
            ia = sort(ia);%need to re-order ia so that things are appropriately listed in order. 
            
            load((waveName))
            fullMaster.(fieldNames{i}).(desigNames{j}).OldWaves = fullMaster.(fieldNames{i}).(desigNames{j}).AllWaves;
            fullMaster.(fieldNames{i}).(desigNames{j}).AllWaves = [];
            fullMaster.(fieldNames{i}).(desigNames{j}).AllWaves = waves(:,:,ia);
            fullMaster.(fieldNames{i}).(desigNames{j}).AverageWaveForms = [];
            fullMaster.(fieldNames{i}).(desigNames{j}).AverageWaveForms = mean(waves(:,:,ia),3);
            fullMaster.(fieldNames{i}).(desigNames{j}).WaveWarning = waveWarn;
            clear 'waves'
        end
    end
end