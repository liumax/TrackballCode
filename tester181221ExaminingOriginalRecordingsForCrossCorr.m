

% clear
targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);

numFiles = length(targetFiles);
rasterWindow = [-0.01 0.01];
rasterVector = [-0.01:0.0005:0.01];

for i = 1:numFiles
    load(targetFiles{i})
    testData = masterData;
    %remove NANs 
%     nanVals = find(isnan(masterData(:,7)));
%     testData(nanVals,:) = [];
    %pull cell types
%     cellTypes = masterHeader(:,7);
    findPVs = find(testData(:,7) == 1);
    findMSNs = find(testData(:,7) == 0);
    probePos = floor(-testData(:,1));
    %if find MSNs and PVs
    if length(findPVs) > 0 & length(findMSNs) > 0
        disp('Both MSNs and PVs found')
        %now we want to analyze for whether MSNs and PVs are close
        %together. 
        %since we want to look within a 400 um radius, we can look across
        %both shanks. 
        disp(strcat(num2str(length(findPVs)),'-FSIs found'))
        for j = 1:length(findPVs)
            distanceVals = sqrt((250*(testData(:,2) - testData(findPVs(j),2))).^2 + (25*(probePos - probePos(findPVs(j)))).^2);
            %now lets pull out the cell types
            msnDistances = distanceVals(findMSNs);
            pvDistances = distanceVals(findPVs);
            %now lets generate a plot. 
            %lets find the nearest square value
            squareVal = ceil(sqrt(size(testData,1)));
            subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.02], [0.02 0.02], [0.01 0.01]);
            hFig = figure;
            set(hFig, 'Position', [10 80 1000 1000])
            format short
            for k = 1:size(testData,1)
                disp(strcat(num2str(k),'/',num2str(size(testData,1))))
                subplot(squareVal,squareVal,k)
                hold on
                alignSpikes = s.(s.DesignationName{findPVs(j)}).SpikeTimes;
                tarSpikes = s.(s.DesignationName{k}).SpikeTimes;
                [rasters] = functionBasicRaster(tarSpikes,alignSpikes,rasterWindow);
%                 [shiftRasters] = functionBasicRaster(tarSpikes+.1,alignSpikes,rasterWindow);
%                 for l = 1:length(alignSpikes)
%                     tmpHist(:,l) = hist(rasters(rasters(:,2) == l,1),rasterVector);
%                 end
                
                histStore = hist(rasters(:,1),rasterVector);
                plot(rasterVector,histStore,'k','LineWidth',2)
                
%                 hist(rasters(:,1),rasterVector)
                plot([0 0],[0 max(histStore)],'r')
                if max(histStore) < 1
                    ylim([0 1])
                else
                    ylim([0 max(histStore)])
                end
                xlim([-0.01 0.01])
                set(gca,'xtick',[])
                if ismember(k,findMSNs)
                    title(strcat(s.DesignationName{k},'-',num2str(round(distanceVals(k)))),'Color','k')
                elseif ismember(k,findPVs) & ~ismember(k,findPVs(j))
                    title(strcat(s.DesignationName{k},'-',num2str(round(distanceVals(k)))),'Color','r')
                elseif ismember(k,findPVs(j))
                    title(strcat(s.DesignationName{k},'TARGET'),'Color','r')
                else
%                     title(strcat(s.DesignationName{k},'-',num2str(distanceVals(k))),'Color','k')
                end
            end
            origName = targetFiles{i};
            spikeGraphName = strcat(origName(1:end-4),s.DesignationName{findPVs(j)},'CrossCorr');
            savefig(hFig,spikeGraphName);

            %save as PDF with correct name
            set(hFig,'Units','Inches');
            pos = get(hFig,'Position');
            set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            print(hFig,spikeGraphName,'-dpdf','-r0')
        end
        %first, determine the number of shanks
        numShanks = length(unique(testData(:,2)));
        for j = 1:length(numShanks)
            disp(strcat('Examining Shank',num2str(j)))
            %lets look for things within a 300 um radius. 
        end
    end
end