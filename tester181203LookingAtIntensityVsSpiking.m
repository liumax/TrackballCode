%lets look at overall plots!  with tone period.
testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
testNames = testNames(findString);

tonePeriod = 2;
numFiles = length(testNames);
counter = 1;
for j = 1:numFiles
    load(testNames{j})
    disp(strcat('Loading-',testNames{j}))
    numUnits = length(s.DesignationName);
    sigThresh = 0.05;
    valThresh = 8;
    sigValStore = [];
    for i = 1:numUnits
        %first, we need to see whether there is a significant response
        SigVals = s.(s.DesignationName{i}).BinSigVals(:,:,tonePeriod);
        %threshold it!
        SigVals(SigVals < sigThresh) = 0;
        SigVals(SigVals >0) = 1;
        SigValsLaser = s.(s.DesignationName{i}).BinSigValsLaser(:,:,tonePeriod);
        %threshold it!
        SigValsLaser(SigValsLaser < sigThresh) = 0;
        SigValsLaser(SigValsLaser >0) = 1;
        numSigVals = length(find(SigVals == 0));
        %since significance doesnt account for positive vs negative
        %findings, lets exclude negative results
        meanBinDiff = mean(s.(s.DesignationName{i}).BinDiff(find(SigVals == 0)));
%         sigValStore(i) = numSigVals;
        %now lets use this as a cutoff
        if numSigVals > valThresh & meanBinDiff > 0
            disp('Significant Crossing!')
            unitName{counter} = strcat(testNames{j},s.DesignationName{i});
            %now lets find the best frequency based on 70db response.
            [Y,I] = max(s.(s.DesignationName{i}).BinTone(:,9)); %I is best frequency
%             [Y,I] = max(mean(s.(s.DesignationName{i}).BinTone')); %I is best frequency
%             mean(s.nt7cluster3.BinTone')
            
            if I
                bfStore(counter) = I;
                intensityResponse(:,i) = s.(s.DesignationName{i}).BinDiff(I,:,tonePeriod);
                intensityResponseLaser(:,i) = s.(s.DesignationName{i}).BinDiffLaser(I,:,tonePeriod);
                bigInt(:,counter) = s.(s.DesignationName{i}).BinDiff(I,:,tonePeriod);
                bigIntLaser(:,counter) = s.(s.DesignationName{i}).BinDiffLaser(I,:,tonePeriod);
                bigSig(:,counter) = SigVals(I,:);
                bigSigLaser(:,counter) = SigValsLaser(I,:);
                fullSig(:,:,counter) = SigVals;
                fullSigTone(:,:,counter) = s.(s.DesignationName{i}).BinSigVals(:,:,tonePeriod);
                fullSigLaser(:,:,counter) = SigValsLaser;
                fullVals(:,:,counter) = s.(s.DesignationName{i}).BinDiff(:,:,tonePeriod);
                fullValsLaser(:,:,counter) = s.(s.DesignationName{i}).BinDiffLaser(:,:,tonePeriod);
                cellType(counter) = masterData(i,7);
                counter = counter + 1;
                tempData = squeeze(s.(s.DesignationName{i}).BinSigVals(:,:,tonePeriod));
                minData = min(tempData);
                minData = find(minData < 0.05,1,'first');
                if minData
                    baseMin(counter) = minData;
                else
                    baseMin(counter) = 0;
                end
                tempData = squeeze(s.(s.DesignationName{i}).BinSigValsLaser(2:end,:,tonePeriod));
                minData = min(tempData);
                minData = find(minData < 0.05,1,'first');
                if minData
                    baseMinLaser(counter) = minData;
                else
                    baseMinLaser(counter) = 0;
                end
%                 figure
%                 subplot(2,1,1)
%                 plot(intensityResponse(:,i),'k')
%                 hold on
%                 plot(intensityResponseLaser(:,i) ,'g')
%                 
%                 subplot(2,1,2)
%                 plot(intensityResponseLaser(:,i) - intensityResponse(:,i))
            end
        end
    end

end

fsis = find(cellType == 1);
msns = find(cellType == 0);

for i = 1:length(bigSig)
    finder = find(bigSig(:,i) == 1,1,'last');
        if finder
        tester(i) = finder;
        end
end

%lets do a basic plot of average response of laser vs no laser. 
for i = 1:length(bigSig)
    %find all significant responses, control
    findSigs = find(fullSig(:,:,i) == 0);
    %pull binned values from this. 
    baseVals = fullVals(:,:,i);
    baseVals = baseVals(findSigs);
    laserVals = fullValsLaser(:,:,i);
    laserVals = laserVals(findSigs);
    signRankResult(i) = signrank(baseVals,laserVals);
    storeVals(i,1) = mean(laserVals);
    storeVals(i,2) = mean(baseVals);
end


hFig = figure;
set(hFig, 'Position', [10 80 600 1000])
subplot(2,1,1)
plot(storeVals(msns,2),storeVals(msns,1),'r.')
hold on
plot([0 max(max(storeVals(msns,:)))],[0 max(max(storeVals(msns,:)))],'k')
title(strcat('Average MSN Responses signrankpval',num2str(signrank(storeVals(msns,2),storeVals(msns,1)))))
xlim([0 max(max(storeVals(msns,:)))])
ylim([0 max(max(storeVals(msns,:)))])
axis square
subplot(2,1,2)
plot(storeVals(fsis,2),storeVals(fsis,1),'r.')

hold on
plot([0 max(max(storeVals(fsis,:)))],[0 max(max(storeVals(fsis,:)))],'k')
title(strcat('Average FSI Responses signrankpval',num2str(signrank(storeVals(fsis,2),storeVals(fsis,1)))))
xlim([0 max(max(storeVals(fsis,:)))])
ylim([0 max(max(storeVals(fsis,:)))])
axis square

spikeGraphName = 'toneBinMSNFSIAverageResponseMagnitudeChange';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')




hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
plot(bigInt(:,fsis))
title('FSI BF Amp Responses')
subplot(2,1,2)
plot(bigInt(:,msns))
title('MSN BF Amp Responses')

spikeGraphName = 'toneBinMSNFSIBFampResponses';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
hold on
plot(mean(bigInt(:,fsis)'),'r')
plot(mean(bigInt(:,msns)'),'k')
plot(mean(bigIntLaser(:,msns)'),'g')
plot(mean(bigIntLaser(:,fsis)'),'m')
title('Mean Amp Resp MSNk FSIr')

subplot(2,1,2)
hold on
plot(mean(bigIntLaser(:,fsis)')-mean(bigInt(:,fsis)'),'r')
plot(mean(bigIntLaser(:,msns)')-mean(bigInt(:,msns)'),'k')
title('Mean Diff Amp Resp MSNk FSIr')

spikeGraphName = 'toneBinMSNFSIPopAvAmpResponses';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
hold on
for i = 1:length(fsis)
    plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('Ind FSI Amp Resp Scatter Laser vs NoLaser')
xlabel('spikes no laser')
ylabel('spikes with laser')
subplot(2,1,2)
hold on
for i = 1:length(msns)
    plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('Ind MSN Amp Resp Scatter Laser vs NoLaser')
xlabel('spikes no laser')
ylabel('spikes with laser')


spikeGraphName = 'toneBinMSNFSIIndivSpikeScatterLaserNoLaser';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
%lets look at average effects for FSIs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStoreFSI = [];
intStoreFSI = [];
valueStore=[];
for i = 1:length(fsis)
    [b,bintr,bintjm] = gmregress(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,fsis(i)),bigIntLaser(:,fsis(i))];
    counter = counter + size(bigInt,1);
    slopeStoreFSI(i) = b(2);
    intStoreFSI(i) = b(1);
end

%determine if there are any values where slope is NaN
nanFind = find(isnan(slopeStoreFSI));
if nanFind
    disp('Some GMREGRESS FAIL FSIs')
    for i = 1:length(nanFind)
        %determine if all laser values are zero
        laserZeroFind = any(bigIntLaser(:,fsis(nanFind(i))));
        %determine if all nonlaser values are zero
        nolaserZeroFind = any(bigInt(:,fsis(nanFind(i))));
        if laserZeroFind == 0 && nolaserZeroFind == 0
            disp('No Spikes, Not adjusting')
        elseif laserZeroFind == 0 && nolaserZeroFind ~= 0
            disp('FULL SILENCING, ADJUSTING SLOPE')
            slopeStoreFSI(nanFind(i)) = 0;
            intStoreFSI(nanFind(i)) = 0;
        end
    end
end

%now we can generate an average based on the slopes and also based on all
%the points

%exclude negative slopes, as these are likely fucked. 
intStoreFSI(slopeStoreFSI<0) = NaN;
slopeStoreFSI(slopeStoreFSI < 0) = NaN;


avSlopeFSI = nanmean(slopeStoreFSI)
avIntFSI = nanmean(intStoreFSI)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
hold on
for i = 1:length(fsis)
    if ~isnan(slopeStoreFSI(i))
        plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.')
        plot(sort(bigInt(:,fsis(i))),sort(bigInt(:,fsis(i)))*slopeStoreFSI(i) + intStoreFSI(i),'k')
    end
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[avIntFSI maxVal*avSlopeFSI],'g')
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('FSI No Laser vs Laser')

spikeGraphName = 'toneBinFSINoLaserVsLaserScatterWithAverageLine';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets look at average effects for MSNs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStore = [];
intStore = [];
valueStore=[];
for i = 1:length(msns)
    [b,bintr,bintjm] = gmregress(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,msns(i)),bigIntLaser(:,msns(i))];
    counter = counter + size(bigInt,1);
    slopeStore(i) = b(2);
    intStore(i) = b(1);
end

%now we can generate an average based on the slopes and also based on all
%the points

intStore(slopeStore<0) = NaN;
slopeStore(slopeStore < 0) = NaN;
intStore(slopeStore>3) = NaN;
slopeStore(slopeStore >3) = NaN;

avSlopeMSN = nanmean(slopeStore)
avIntMSN = nanmean(intStore)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);


hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
hold on
for i = 1:length(msns)
    if ~isnan(slopeStore(i))
        plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.')
        plot(bigInt(:,msns(i)),bigInt(:,msns(i))*slopeStore(i) + intStore(i),'k')
    end
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[avIntMSN maxVal*avSlopeMSN],'g')
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('MSN No Laser vs Laser')


spikeGraphName = 'toneBinMSNNoLaserVsLaserScatterWithAverageLine';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets plot out distribution of slopes. 
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,2,1)
hold on
hist(slopeStore,[0:0.1:3])
tester = slopeStore;
tester(isnan(tester)) = [];
plot(nanmedian(tester),1,'g.')
plot(nanmean(tester),1.5,'r.')
title(strcat('MSN Slope Dist median green mean red mean',num2str(nanmean(slopeStore)), 'pval',num2str(signrank(slopeStore-1))))
xlim([0 3])
subplot(2,2,2)
hold on
hist(intStore,20)
plot(nanmedian(intStore),1,'g.')
plot(nanmean(intStore),1.5,'r.')
title(strcat('MSN Int Dist median green mean red mean',num2str(nanmean(intStore)), ' pval',num2str(signrank(intStore))))
subplot(2,2,3)
hold on
hist(slopeStoreFSI,[0:0.1:3])
plot(nanmedian(slopeStoreFSI),1,'g.')
plot(nanmean(slopeStoreFSI),1.5,'r.')
title(strcat('FSI Slope Dist median green mean red mean',num2str(nanmean(slopeStoreFSI)), ' pval',num2str(signrank(slopeStoreFSI-1))))
xlim([0 3])
subplot(2,2,4)
hold on
hist(intStoreFSI,20)
plot(nanmedian(intStoreFSI),1,'g.')
plot(nanmean(intStoreFSI),1.5,'r.')
title(strcat('FSI Int Dist median green mean red mean',num2str(nanmean(intStoreFSI)), ' pval',num2str(signrank(intStoreFSI))))

spikeGraphName = 'toneBinMSNFSIDistributionSlopesAndIntercepts';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets try alternative where we average across multiple stimuli

%lets do a basic plot of average response of laser vs no laser. 
for i = 1:length(bigSig)
    %find all significant responses, control
    findSigs = find(fullSig(:,:,i) == 0);
    %pull binned values from this. 
    baseVals = fullVals(:,:,i);
    baseVals = baseVals(findSigs);
    laserVals = fullValsLaser(:,:,i);
    laserVals = laserVals(findSigs);
    [b,bintr,bintjm] = gmregress(baseVals,laserVals,sigVal);
    combSlopeStore(i) = b(2);
    combIntStore(i) = b(1);
end

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,2,1)
hold on
hist(combSlopeStore(msns),[0:0.1:3])
tester = combSlopeStore(msns);
tester(isnan(tester)) = [];
plot(nanmedian(tester),1,'g.')
plot(nanmean(tester),1.5,'r.')
title(strcat('MSN Slope Dist median green mean red mean',num2str(nanmean(combSlopeStore(msns))), 'pval',num2str(signrank(combSlopeStore(msns)-1))))
xlim([0 3])
subplot(2,2,2)
hold on
hist(combIntStore(msns),20)
plot(nanmedian(combIntStore(msns)),1,'g.')
plot(nanmean(combIntStore(msns)),1.5,'r.')
title(strcat('MSN Int Dist median green mean red mean',num2str(nanmean(combIntStore(msns))), ' pval',num2str(signrank(combIntStore(msns)))))
subplot(2,2,3)
hold on
hist(combSlopeStore(fsis),[0:0.1:3])
plot(nanmedian(combSlopeStore(fsis)),1,'g.')
plot(nanmean(combSlopeStore(fsis)),1.5,'r.')
title(strcat('FSI Slope Dist median green mean red mean',num2str(nanmean(combSlopeStore(fsis))), ' pval',num2str(signrank(combSlopeStore(fsis)-1))))
xlim([0 3])
subplot(2,2,4)
hold on
hist(combIntStore(fsis),20)
plot(nanmedian(combIntStore(fsis)),1,'g.')
plot(nanmean(combIntStore(fsis)),1.5,'r.')
title(strcat('FSI Int Dist median green mean red mean',num2str(nanmean(combIntStore(fsis))), ' pval',num2str(signrank(combIntStore(fsis)))))

spikeGraphName = 'toneBinMSNFSIAllSigRespMultiFreqDistributionSlopesAndIntercepts';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets try and do this thing where we find thresholds etc. First thing
%we need to do with this is clean up the data. 
fullSigTone(fullSigTone < 0.05) = 0;
fullSigTone(fullSigTone >0) = 1;
% for i =1:length(msns)
%     figure
%     subplot(2,1,1)
%     imagesc(fullSigTone(:,:,msns(i)))
%     colormap('parula')
%     title(unitName{i}, 'Interpreter', 'none')
%     subplot(2,1,2)
%     imagesc(fullSigLaser(:,:,msns(i)))
%     colormap('parula')
% end

%lets try and clean up this data, by eliminating things that arent
%connected to highest amplitude. This means going from the end, then
%marching down until I hit the first non-significant value. Remember that
%0s are significant, 1s are not

fullSigToneAnnot = fullSigTone;
fullSigLaserAnnot = fullSigLaser;

annotDBVals = zeros(length(fullSig),size(fullSigTone,1));
annotDBValsLaser = zeros(length(fullSig),size(fullSigTone,1));

for i = 1:length(fullSig)
    testData = fullSigTone(:,:,i);
    for j = 1:size(testData,1) %going by freq
        whileTrig = 0;
        whileCount = size(testData,2);
        while whileTrig == 0
            %determine if first value is zero or not
            if testData(j,whileCount) == 0
                disp(strcat('Significant Response at DB',num2str(whileCount)))
                if whileCount == 1
                    disp('Hit Lowest, Moving On to Next Freq')
                    annotDBVals(i,j) = whileCount;
                    whileTrig = 1;
                else
                    disp('Moving to Lower Amp')
                    whileCount = whileCount -1;
                end
            elseif testData(j,whileCount) == 1
                disp(strcat('First Non-Sig Resp at DB',num2str(whileCount)))
                annotDBVals(i,j) = whileCount+1;
                testData(j,1:whileCount) = 1;
                whileTrig = 1;
            end
        end
    end
    fullSigToneAnnot(:,:,i) = testData;
    disp('Finished a unit!')
end

for i = 1:length(fullSig)
    testData = fullSigLaser(:,:,i);
    for j = 1:size(testData,1) %going by freq
        whileTrig = 0;
        whileCount = size(testData,2);
        while whileTrig == 0
            %determine if first value is zero or not
            if testData(j,whileCount) == 0
                disp(strcat('Significant Response at DB',num2str(whileCount)))
                if whileCount == 1
                    disp('Hit Lowest, Moving On to Next Freq')
                    annotDBValsLaser(i,j) = whileCount;
                    whileTrig = 1;
                else
                    disp('Moving to Lower Amp')
                    whileCount = whileCount -1;
                end
            elseif testData(j,whileCount) == 1
                disp(strcat('First Non-Sig Resp at DB',num2str(whileCount)))
                annotDBValsLaser(i,j) = whileCount+1;
                testData(j,1:whileCount) = 1;
                whileTrig = 1;
            end
        end
    end
    fullSigLaserAnnot(:,:,i) = testData;
    disp('Finished a unit!')
end

% for i =1:length(msns)
%     figure
%     subplot(2,1,1)
%     imagesc(fullSigToneAnnot(:,:,msns(i)))
%     colormap('parula')
%     title(unitName{msns(i)}, 'Interpreter', 'none')
%     subplot(2,1,2)
%     imagesc(fullSigLaserAnnot(:,:,msns(i)))
%     colormap('parula')
% end

%this seems to work! Now lets eliminate units without enough responses
%anymore

for i = 1:length(fullSig)
    AnnotVals(i) = sum(sum(fullSigToneAnnot(:,:,i)));
end

figure
hist(AnnotVals,100)

%based on this, lets eliminate things with values greater than or equal to
%50

findBadAnnot = find(AnnotVals >=50);

annotMSNs = setdiff(msns,findBadAnnot);
annotFSIs = setdiff(fsis,findBadAnnot);

annotMSNDiffVals = annotDBVals(annotMSNs,:) -  annotDBValsLaser(annotMSNs,:);
annotFSIDiffVals = annotDBVals(annotFSIs,:) -  annotDBValsLaser(annotFSIs,:);

%calculate mean changes per unit
meanAnnotMSNDiff = mean(annotMSNDiffVals');
meanAnnotFSIDiff = mean(annotFSIDiffVals');

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
hist(meanAnnotMSNDiff,10)
title(strcat('MSN Mean Difference in DB Threshold mean',num2str(mean(meanAnnotMSNDiff)), 'pvalsignrank',num2str(signrank(mean(annotDBVals(annotMSNs,:)'),mean(annotDBValsLaser(annotMSNs,:)')))))
subplot(2,1,2)
hist(meanAnnotFSIDiff,10)
title(strcat('FSI Mean Difference in DB Threshold mean',num2str(mean(meanAnnotFSIDiff)), ' pvalsignrank',num2str(signrank(mean(annotDBVals(annotFSIs,:)'),mean(annotDBValsLaser(annotFSIs,:)')))))

spikeGraphName = 'toneBinFSIandMSNMeanDBThreshChange';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets try and do the same thing based on just the first threshold
%response per unit. 
fullSigTone;
fullSigLaser;
for i = 1:length(fullSigTone)
    for j = 1:size(fullSigTone,1)
        finder = find(fullSigTone(j,:,i) == 0,1,'first');
        if finder
            absDBthreshVal(i,j) = finder;
        end
        finder = find(fullSigLaser(j,:,i) == 0,1,'first');
        if finder
            absDBthreshValLaser(i,j) = finder;
        end
    end
end

%find MSNs
meanAbsDBThreshChange = mean(absDBthreshVal') - mean(absDBthreshValLaser');

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
hist(meanAbsDBThreshChange(msns),10)
title(strcat('MSN Mean Difference in Absolute DB Threshold mean',num2str(mean(meanAbsDBThreshChange(msns))), 'pvalsignrank',num2str(signrank(mean(absDBthreshVal(msns,:)'),mean(absDBthreshValLaser(msns,:)')))))
subplot(2,1,2)
hist(meanAbsDBThreshChange(fsis),10)
title(strcat('FSI Mean Difference in Absolute DB Threshold mean',num2str(mean(meanAbsDBThreshChange(fsis))), ' pvalsignrank',num2str(signrank(mean(absDBthreshVal(fsis,:)'),mean(absDBthreshValLaser(fsis,:)')))))

spikeGraphName = 'toneBinFSIandMSNABSOLUTEMeanDBThreshChange';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%Now this method isnt exactly great either: predisposed to noise etc. 


%lets try an alternative. We could use BF, and if BF doesnt work, we can
%shift left or right. 
fullSigTone;
fullSigLaser;
for i = 1:length(bfStore)
    disp('NewUnit!')
    whileTrig = 0;
%     whileTrig2 = 0
    tarFreqs = [1:1:size(fullSigTone,1)];
    tarFreq = bfStore(i);
    while whileTrig == 0;
        if length(tarFreqs) == 0
            selThresh(i,1) = NaN;
            selThresh(i,2) = NaN;
            break
        end
        %pull the BF significance
        bfLine = fullSigTone(tarFreq,:,i);
        %if the bfLine is just zeros, or threshold value is beyond the first db
        %step, then we need to move left or right.
        if sum(bfLine) == size(fullSigTone,2)
            disp('FAILURE: NO SIGNIFICANT VALUES AT BF')
%             while whileTrig2 == 0
            tarFreqs(tarFreqs == tarFreq) = [];
            %find difference between BF and current freq list
            diffFind = abs(tarFreqs - bfStore(i));
            %find closest to original BF
            [B I] = min(diffFind);
            %set this as the new tarFreq
            tarFreq = tarFreqs(I);
%             end
        elseif bfLine(1) == 0 && sum(bfLine) < size(fullSigTone,2)
            disp('Best Freq Has Maxed Threshold')
            if bfStore(i) == 1
                
                tarFreqs(tarFreqs == tarFreq) = [];
                tarFreq = tarFreqs(1);
            else
                tarFreqs(tarFreqs == tarFreq) = [];
                %find difference between BF and current freq list
                diffFind = abs(tarFreqs - bfStore(i));
                %find closest to original BF
                [B I] = min(diffFind);
                %set this as the new tarFreq
                tarFreq = tarFreqs(I);
            end
        elseif bfLine(1) ~= 0 && sum(bfLine) < size(fullSigTone,2)
            %if best freq doesnt have threshold at first value, then save and
            %record threshold value for laser trials
            disp('Best Freq Thresh Found')
            selThresh(i,1) = find(bfLine == 0,1,'first');
            if find(fullSigLaser(tarFreq,:,i) == 0,1,'first')
                selThresh(i,2) = find(fullSigLaser(tarFreq,:,i) == 0,1,'first');
            else
                selThresh(i,2) = NaN;
            end
            break
        end
    end
end

%now lets pull msns vs fsis.
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
hist(selThresh(msns,1)-selThresh(msns,2),10)
title(strcat('MSN Mean Difference in BF/Alt DB Threshold mean',num2str(nanmean(selThresh(msns,1)-selThresh(msns,2))), 'pvalsignrank',num2str(signrank(selThresh(msns,1),selThresh(msns,2)))))
subplot(2,1,2)
hist(selThresh(fsis,1)-selThresh(fsis,2),10)
title(strcat('FSI Mean Difference in BF/Alt DB Threshold mean',num2str(nanmean(selThresh(fsis,1)-selThresh(fsis,2))), ' pvalsignrank',num2str(signrank(selThresh(fsis,1),selThresh(fsis,2)))))

spikeGraphName = 'toneBinFSIandMSNBFMeanDBThreshChange';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now for a sanity check: lets see what the slopes are for the BF binned
%datasets. Possibly, we could select for ones that are largely increasing?
for i = 1:length(bigInt)
    p = polyfit([1:1:9]',bigInt(:,i),1);
    fitVal(i,:) = p;
end

bigFits = find(fitVal > 0.1);


% signrank(meanAnnotFSIDiff)

%% I"m a buffer!

clear
%% Lets try without baseline subtraction
%lets look at overall plots!  with tone period.
testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
testNames = testNames(findString);

tonePeriod = 2;
numFiles = length(testNames);
counter = 1;
for j = 1:numFiles
    load(testNames{j})
    disp(strcat('Loading-',testNames{j}))
    numUnits = length(s.DesignationName);
    sigThresh = 0.05;
    valThresh = 8;
    sigValStore = [];
    for i = 1:numUnits
        %first, we need to see whether there is a significant response
        SigVals = s.(s.DesignationName{i}).BinSigVals(:,:,tonePeriod);
        %threshold it!
        SigVals(SigVals < sigThresh) = 0;
        SigVals(SigVals >0) = 1;
        SigValsLaser = s.(s.DesignationName{i}).BinSigValsLaser(:,:,tonePeriod);
        %threshold it!
        SigValsLaser(SigValsLaser < sigThresh) = 0;
        SigValsLaser(SigValsLaser >0) = 1;
        numSigVals = length(find(SigVals == 0));
        %since significance doesnt account for positive vs negative
        %findings, lets exclude negative results
        meanBinDiff = mean(s.(s.DesignationName{i}).BinDiff(find(SigVals == 0)));
%         sigValStore(i) = numSigVals;
        %now lets use this as a cutoff
        if numSigVals > valThresh & meanBinDiff > 0
            disp('Significant Crossing!')
            unitName{counter} = strcat(testNames{j},s.DesignationName{i});
            %now lets find the best frequency based on 70db response.
            [Y,I] = max(s.(s.DesignationName{i}).BinTone(:,9)); %I is best frequency
%             [Y,I] = max(mean(s.(s.DesignationName{i}).BinTone')); %I is best frequency
%             mean(s.nt7cluster3.BinTone')
            
            if I
                bfStore(counter) = I;
                intensityResponse(:,i) = s.(s.DesignationName{i}).BinTone(I,:);
                intensityResponseLaser(:,i) = s.(s.DesignationName{i}).BinToneLaser(I,:);
                bigInt(:,counter) = s.(s.DesignationName{i}).BinTone(I,:);
                bigIntLaser(:,counter) = s.(s.DesignationName{i}).BinToneLaser(I,:);
                bigSig(:,counter) = SigVals(I,:);
                bigSigLaser(:,counter) = SigValsLaser(I,:);
                fullSig(:,:,counter) = SigVals;
                fullSigTone(:,:,counter) = s.(s.DesignationName{i}).BinSigVals(:,:,2);
                fullSigLaser(:,:,counter) = SigValsLaser;
                fullVals(:,:,counter) = s.(s.DesignationName{i}).BinTone(:,:);
                fullValsLaser(:,:,counter) = s.(s.DesignationName{i}).BinToneLaser(:,:);
                cellType(counter) = masterData(i,7);
                counter = counter + 1;
                tempData = squeeze(s.(s.DesignationName{i}).BinSigVals(:,:,tonePeriod));
                minData = min(tempData);
                minData = find(minData < 0.05,1,'first');
                if minData
                    baseMin(counter) = minData;
                else
                    baseMin(counter) = 0;
                end
                tempData = squeeze(s.(s.DesignationName{i}).BinSigValsLaser(2:end,:,tonePeriod));
                minData = min(tempData);
                minData = find(minData < 0.05,1,'first');
                if minData
                    baseMinLaser(counter) = minData;
                else
                    baseMinLaser(counter) = 0;
                end
%                 figure
%                 subplot(2,1,1)
%                 plot(intensityResponse(:,i),'k')
%                 hold on
%                 plot(intensityResponseLaser(:,i) ,'g')
%                 
%                 subplot(2,1,2)
%                 plot(intensityResponseLaser(:,i) - intensityResponse(:,i))
            end
        end
    end

end

fsis = find(cellType == 1);
msns = find(cellType == 0);

for i = 1:length(bigSig)
    finder = find(bigSig(:,i) == 1,1,'last');
        if finder
        tester(i) = finder;
        end
end

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
plot(bigInt(:,fsis))
title('FSI BF Amp Responses')
subplot(2,1,2)
plot(bigInt(:,msns))
title('MSN BF Amp Responses')

spikeGraphName = 'NoSUBMSNFSIBFampResponses';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
hold on
plot(mean(bigInt(:,fsis)'),'r')
plot(mean(bigInt(:,msns)'),'k')
plot(mean(bigIntLaser(:,msns)'),'g')
plot(mean(bigIntLaser(:,fsis)'),'m')
title('Mean Amp Resp MSNk FSIr')

subplot(2,1,2)
hold on
plot(mean(bigIntLaser(:,fsis)')-mean(bigInt(:,fsis)'),'r')
plot(mean(bigIntLaser(:,msns)')-mean(bigInt(:,msns)'),'k')
title('Mean Diff Amp Resp MSNk FSIr')

spikeGraphName = 'NoSUBMSNFSIPopAvAmpResponses';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
hold on
for i = 1:length(fsis)
    plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('Ind FSI Amp Resp Scatter Laser vs NoLaser')
xlabel('spikes no laser')
ylabel('spikes with laser')
subplot(2,1,2)
hold on
for i = 1:length(msns)
    plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('Ind MSN Amp Resp Scatter Laser vs NoLaser')
xlabel('spikes no laser')
ylabel('spikes with laser')


spikeGraphName = 'NoSUBMSNFSIIndivSpikeScatterLaserNoLaser';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
%lets look at average effects for FSIs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStoreFSI = [];
intStoreFSI = [];
valueStore=[];
for i = 1:length(fsis)
    [b,bintr,bintjm] = gmregress(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,fsis(i)),bigIntLaser(:,fsis(i))];
    counter = counter + size(bigInt,1);
    slopeStoreFSI(i) = b(2);
    intStoreFSI(i) = b(1);
end

%determine if there are any values where slope is NaN
nanFind = find(isnan(slopeStoreFSI));
if nanFind
    disp('Some GMREGRESS FAIL FSIs')
    for i = 1:length(nanFind)
        %determine if all laser values are zero
        laserZeroFind = any(bigIntLaser(:,fsis(nanFind(i))));
        %determine if all nonlaser values are zero
        nolaserZeroFind = any(bigInt(:,fsis(nanFind(i))));
        if laserZeroFind == 0 && nolaserZeroFind == 0
            disp('No Spikes, Not adjusting')
        elseif laserZeroFind == 0 && nolaserZeroFind ~= 0
            disp('FULL SILENCING, ADJUSTING SLOPE')
            slopeStoreFSI(nanFind(i)) = 0;
            intStoreFSI(nanFind(i)) = 0;
        end
    end
end

%now we can generate an average based on the slopes and also based on all
%the points

%exclude negative slopes, as these are likely fucked. 
intStoreFSI(slopeStoreFSI<0) = NaN;
slopeStoreFSI(slopeStoreFSI < 0) = NaN;


avSlopeFSI = nanmean(slopeStoreFSI)
avIntFSI = nanmean(intStoreFSI)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
hold on
for i = 1:length(fsis)
    if ~isnan(slopeStoreFSI(i))
        plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.')
        plot(sort(bigInt(:,fsis(i))),sort(bigInt(:,fsis(i)))*slopeStoreFSI(i) + intStoreFSI(i),'k')
    end
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[avIntFSI maxVal*avSlopeFSI],'g')
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('FSI No Laser vs Laser')

spikeGraphName = 'NoSUBFSINoLaserVsLaserScatterWithAverageLine';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets look at average effects for MSNs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStore = [];
intStore = [];
valueStore=[];
for i = 1:length(msns)
    [b,bintr,bintjm] = gmregress(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,msns(i)),bigIntLaser(:,msns(i))];
    counter = counter + size(bigInt,1);
    slopeStore(i) = b(2);
    intStore(i) = b(1);
end

%now we can generate an average based on the slopes and also based on all
%the points

intStore(slopeStore<0) = NaN;
slopeStore(slopeStore < 0) = NaN;
intStore(slopeStore>3) = NaN;
slopeStore(slopeStore >3) = NaN;

avSlopeMSN = nanmean(slopeStore)
avIntMSN = nanmean(intStore)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);


hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
hold on
for i = 1:length(msns)
    if ~isnan(slopeStore(i))
        plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.')
        plot(bigInt(:,msns(i)),bigInt(:,msns(i))*slopeStore(i) + intStore(i),'k')
    end
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[avIntMSN maxVal*avSlopeMSN],'g')
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('MSN No Laser vs Laser')


spikeGraphName = 'NoSUBMSNNoLaserVsLaserScatterWithAverageLine';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets plot out distribution of slopes. 
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,2,1)
hold on
hist(slopeStore,20)
tester = slopeStore;
tester(isnan(tester)) = [];
plot(nanmedian(tester),1,'g.')
plot(nanmean(tester),1.5,'r.')
title(strcat('MSN Slope Dist median green mean red mean',num2str(nanmean(slopeStore)), 'pval',num2str(signrank(slopeStore))))
subplot(2,2,2)
hold on
hist(intStore,20)
plot(nanmedian(intStore),1,'g.')
plot(nanmean(intStore),1.5,'r.')
title(strcat('MSN Int Dist median green mean red mean',num2str(nanmean(intStore)), ' pval',num2str(signrank(intStore))))
subplot(2,2,3)
hold on
hist(slopeStoreFSI,20)
plot(nanmedian(slopeStoreFSI),1,'g.')
plot(nanmean(slopeStoreFSI),1.5,'r.')
title(strcat('FSI Slope Dist median green mean red mean',num2str(nanmean(slopeStoreFSI)), ' pval',num2str(signrank(slopeStoreFSI))))
subplot(2,2,4)
hold on
hist(intStoreFSI,20)
plot(nanmedian(intStoreFSI),1,'g.')
plot(nanmean(intStoreFSI),1.5,'r.')
title(strcat('FSI Int Dist median green mean red mean',num2str(nanmean(intStoreFSI)), ' pval',num2str(signrank(intStoreFSI))))

spikeGraphName = 'NoSUBMSNFSIDistributionSlopesAndIntercepts';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets try and do this thing where we find thresholds etc. First thing
%we need to do with this is clean up the data. 
fullSigTone(fullSigTone < 0.05) = 0;
fullSigTone(fullSigTone >0) = 1;
for i =1:length(msns)
    figure
    subplot(2,1,1)
    imagesc(fullSigTone(:,:,msns(i)))
    colormap('parula')
    title(unitName{i}, 'Interpreter', 'none')
    subplot(2,1,2)
    imagesc(fullSigLaser(:,:,msns(i)))
    colormap('parula')
end

%lets try and clean up this data, by eliminating things that arent
%connected to highest amplitude. This means going from the end, then
%marching down until I hit the first non-significant value. Remember that
%0s are significant, 1s are not

fullSigToneAnnot = fullSigTone;
fullSigLaserAnnot = fullSigLaser;

annotDBVals = zeros(length(fullSig),size(fullSigTone,1));
annotDBValsLaser = zeros(length(fullSig),size(fullSigTone,1));

for i = 1:length(fullSig)
    testData = fullSigTone(:,:,i);
    for j = 1:size(testData,1) %going by freq
        whileTrig = 0;
        whileCount = size(testData,2);
        while whileTrig == 0
            %determine if first value is zero or not
            if testData(j,whileCount) == 0
                disp(strcat('Significant Response at DB',num2str(whileCount)))
                if whileCount == 1
                    disp('Hit Lowest, Moving On to Next Freq')
                    annotDBVals(i,j) = whileCount;
                    whileTrig = 1;
                else
                    disp('Moving to Lower Amp')
                    whileCount = whileCount -1;
                end
            elseif testData(j,whileCount) == 1
                disp(strcat('First Non-Sig Resp at DB',num2str(whileCount)))
                annotDBVals(i,j) = whileCount+1;
                testData(j,1:whileCount) = 1;
                whileTrig = 1;
            end
        end
    end
    fullSigToneAnnot(:,:,i) = testData;
    disp('Finished a unit!')
end

for i = 1:length(fullSig)
    testData = fullSigLaser(:,:,i);
    for j = 1:size(testData,1) %going by freq
        whileTrig = 0;
        whileCount = size(testData,2);
        while whileTrig == 0
            %determine if first value is zero or not
            if testData(j,whileCount) == 0
                disp(strcat('Significant Response at DB',num2str(whileCount)))
                if whileCount == 1
                    disp('Hit Lowest, Moving On to Next Freq')
                    annotDBValsLaser(i,j) = whileCount;
                    whileTrig = 1;
                else
                    disp('Moving to Lower Amp')
                    whileCount = whileCount -1;
                end
            elseif testData(j,whileCount) == 1
                disp(strcat('First Non-Sig Resp at DB',num2str(whileCount)))
                annotDBValsLaser(i,j) = whileCount+1;
                testData(j,1:whileCount) = 1;
                whileTrig = 1;
            end
        end
    end
    fullSigLaserAnnot(:,:,i) = testData;
    disp('Finished a unit!')
end

% for i =1:length(msns)
%     figure
%     subplot(2,1,1)
%     imagesc(fullSigToneAnnot(:,:,msns(i)))
%     colormap('parula')
%     title(unitName{msns(i)}, 'Interpreter', 'none')
%     subplot(2,1,2)
%     imagesc(fullSigLaserAnnot(:,:,msns(i)))
%     colormap('parula')
% end

%this seems to work! Now lets eliminate units without enough responses
%anymore

for i = 1:length(fullSig)
    AnnotVals(i) = sum(sum(fullSigToneAnnot(:,:,i)));
end

figure
hist(AnnotVals,100)

%based on this, lets eliminate things with values greater than or equal to
%50

findBadAnnot = find(AnnotVals >=50);

annotMSNs = setdiff(msns,findBadAnnot);
annotFSIs = setdiff(fsis,findBadAnnot);

annotMSNDiffVals = annotDBVals(annotMSNs,:) -  annotDBValsLaser(annotMSNs,:);
annotFSIDiffVals = annotDBVals(annotFSIs,:) -  annotDBValsLaser(annotFSIs,:);

%calculate mean changes per unit
meanAnnotMSNDiff = mean(annotMSNDiffVals');
meanAnnotFSIDiff = mean(annotFSIDiffVals');

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
hist(meanAnnotMSNDiff,10)
title(strcat('MSN Mean Difference in DB Threshold mean',num2str(mean(meanAnnotMSNDiff)), 'pval',num2str(signrank(meanAnnotMSNDiff))))
subplot(2,1,2)
hist(meanAnnotFSIDiff,10)
title(strcat('FSI Mean Difference in DB Threshold mean',num2str(mean(meanAnnotFSIDiff)), ' pval',num2str(signrank(meanAnnotFSIDiff))))

spikeGraphName = 'NoSUBFSIandMSNMeanDBThreshChange';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets try and do the same thing based on just the first threshold
%response per unit. 
fullSigTone;
fullSigLaser;
for i = 1:length(fullSigTone)
    for j = 1:size(fullSigTone,1)
        finder = find(fullSigTone(j,:,i) == 0,1,'first');
        if finder
            absDBthreshVal(i,j) = finder;
        end
        finder = find(fullSigLaser(j,:,i) == 0,1,'first');
        if finder
            absDBthreshValLaser(i,j) = finder;
        end
    end
end

%find MSNs
meanAbsDBThreshChange = mean(absDBthreshVal') - mean(absDBthreshValLaser');

hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
hist(meanAbsDBThreshChange(msns),10)
title(strcat('MSN Mean Difference in Absolute DB Threshold mean',num2str(mean(meanAbsDBThreshChange(msns))), 'pval',num2str(signrank(meanAbsDBThreshChange(msns)))))
subplot(2,1,2)
hist(meanAbsDBThreshChange(fsis),10)
title(strcat('FSI Mean Difference in Absolute DB Threshold mean',num2str(mean(meanAbsDBThreshChange(fsis))), ' pval',num2str(signrank(meanAbsDBThreshChange(fsis)))))

spikeGraphName = 'NoSUBFSIandMSNABSOLUTEMeanDBThreshChange';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%Now this method isnt exactly great either: predisposed to noise etc. 


%lets try an alternative. We could use BF, and if BF doesnt work, we can
%shift left or right. 
fullSigTone;
fullSigLaser;
for i = 1:length(bfStore)
    disp('NewUnit!')
    whileTrig = 0;
%     whileTrig2 = 0
    tarFreqs = [1:1:size(fullSigTone,1)];
    tarFreq = bfStore(i);
    while whileTrig == 0;
        if length(tarFreqs) == 0
            selThresh(i,1) = NaN;
            selThresh(i,2) = NaN;
            break
        end
        %pull the BF significance
        bfLine = fullSigTone(tarFreq,:,i);
        %if the bfLine is just zeros, or threshold value is beyond the first db
        %step, then we need to move left or right.
        if sum(bfLine) == size(fullSigTone,2)
            disp('FAILURE: NO SIGNIFICANT VALUES AT BF')
%             while whileTrig2 == 0
            tarFreqs(tarFreqs == tarFreq) = [];
            %find difference between BF and current freq list
            diffFind = abs(tarFreqs - bfStore(i));
            %find closest to original BF
            [B I] = min(diffFind);
            %set this as the new tarFreq
            tarFreq = tarFreqs(I);
%             end
        elseif bfLine(1) == 0 && sum(bfLine) < size(fullSigTone,2)
            disp('Best Freq Has Maxed Threshold')
            if bfStore(i) == 1
                
                tarFreqs(tarFreqs == tarFreq) = [];
                tarFreq = tarFreqs(1);
            else
                tarFreqs(tarFreqs == tarFreq) = [];
                %find difference between BF and current freq list
                diffFind = abs(tarFreqs - bfStore(i));
                %find closest to original BF
                [B I] = min(diffFind);
                %set this as the new tarFreq
                tarFreq = tarFreqs(I);
            end
        elseif bfLine(1) ~= 0 && sum(bfLine) < size(fullSigTone,2)
            %if best freq doesnt have threshold at first value, then save and
            %record threshold value for laser trials
            disp('Best Freq Thresh Found')
            selThresh(i,1) = find(bfLine == 0,1,'first');
            if find(fullSigLaser(tarFreq,:,i) == 0,1,'first')
                selThresh(i,2) = find(fullSigLaser(tarFreq,:,i) == 0,1,'first');
            else
                selThresh(i,2) = NaN;
            end
            break
        end
    end
end

%now lets pull msns vs fsis.
hFig = figure;
set(hFig, 'Position', [10 80 1900 1000])
subplot(2,1,1)
hist(selThresh(msns,1)-selThresh(msns,2),10)
title(strcat('MSN Mean Difference in BF/Alt DB Threshold mean',num2str(nanmean(selThresh(msns,1)-selThresh(msns,2))), 'pval',num2str(signrank(selThresh(msns,1)-selThresh(msns,2)))))
subplot(2,1,2)
hist(selThresh(fsis,1)-selThresh(fsis,2),10)
title(strcat('FSI Mean Difference in BF/Alt DB Threshold mean',num2str(nanmean(selThresh(fsis,1)-selThresh(fsis,2))), ' pval',num2str(signrank(selThresh(fsis,1)-selThresh(fsis,2)))))

spikeGraphName = 'NoSUBFSIandMSNBFMeanDBThreshChange';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
% signrank(meanAnnotFSIDiff)

clear

%% now lets try the same with peak values
testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
testNames = testNames(findString);


numFiles = length(testNames);
counter = 1;
for j = 1:numFiles
    load(testNames{j})
    disp(strcat('Loading-',testNames{j}))
    numUnits = length(s.DesignationName);
    sigThresh = 0.05;
    valThresh = 8;
    sigValStore = [];
    for i = 1:numUnits
        %first, we need to see whether there is a significant response
        SigVals = s.(s.DesignationName{i}).BinSigVals(:,:,2);
        %threshold it!
        SigVals(SigVals < sigThresh) = 0;
        SigVals(SigVals >0) = 1;
        numSigVals = length(find(SigVals == 0));
        sigValStore(i) = numSigVals;
        %now lets use this as a cutoff
        if numSigVals > valThresh
            disp('Significant Crossing!')
            unitName{counter} = strcat(testNames{j},s.DesignationName{i});
            %now lets find the best frequency based on 70db response.
            [Y,I] = max(s.(s.DesignationName{i}).PeakMapTone(:,end)); %I is best frequency
            if I
                intensityResponse(:,i) = s.(s.DesignationName{i}).PeakMapTone(I,:);
                intensityResponseLaser(:,i) = s.(s.DesignationName{i}).PeakMapToneLaser(I,:);
                bigInt(:,counter) = s.(s.DesignationName{i}).PeakMapTone(I,:);
                bigIntLaser(:,counter) = s.(s.DesignationName{i}).PeakMapToneLaser(I,:);
                cellType(counter) = masterData(i,7);
                counter = counter + 1;
%                 figure
%                 subplot(2,1,1)
%                 plot(intensityResponse(:,i),'k')
%                 hold on
%                 plot(intensityResponseLaser(:,i) ,'g')
%                 
%                 subplot(2,1,2)
%                 plot(intensityResponseLaser(:,i) - intensityResponse(:,i))
            end
        end
    end

end

fsis = find(cellType == 1);
msns = find(cellType == 0);

hFig = figure;
subplot(2,1,1)
plot(bigInt(:,fsis))
title('FSI BF Amp Responses')
subplot(2,1,2)
plot(bigInt(:,msns))
title('MSN BF Amp Responses')

spikeGraphName = 'peakBasedMSNFSIBFampResponses';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
subplot(2,1,1)
hold on
plot(mean(bigInt(:,fsis)'),'r')
plot(mean(bigInt(:,msns)'),'k')
plot(mean(bigIntLaser(:,msns)'),'g')
plot(mean(bigIntLaser(:,fsis)'),'m')
title('Mean Amp Resp MSNk FSIr')

subplot(2,1,2)
hold on
plot(mean(bigIntLaser(:,fsis)')-mean(bigInt(:,fsis)'),'r')
plot(mean(bigIntLaser(:,msns)')-mean(bigInt(:,msns)'),'k')
title('Mean Diff Amp Resp MSNk FSIr')

spikeGraphName = 'peakBasedMSNFSIPopAvAmpResponses';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
subplot(2,1,1)
hold on
for i = 1:length(fsis)
    plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('Ind FSI Amp Resp Scatter Laser vs NoLaser')
xlabel('spikes no laser')
ylabel('spikes with laser')

figure
hold on
for i = 1:length(msns)
    plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('Ind MSN Amp Resp Scatter Laser vs NoLaser')
xlabel('spikes no laser')
ylabel('spikes with laser')


spikeGraphName = 'peakBasedMSNFSIIndivSpikeScatterLaserNoLaser';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
%lets look at average effects for FSIs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStoreFSI = [];
intStoreFSI = [];
valueStore=[];
for i = 1:length(fsis)
    [b,bintr,bintjm] = gmregress(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,fsis(i)),bigIntLaser(:,fsis(i))];
    counter = counter + size(bigInt,1);
    slopeStoreFSI(i) = b(2);
    intStoreFSI(i) = b(1);
end

%determine if there are any values where slope is NaN
nanFind = find(isnan(slopeStoreFSI));
if nanFind
    disp('Some GMREGRESS FAIL FSIs')
    for i = 1:length(nanFind)
        %determine if all laser values are zero
        laserZeroFind = any(bigIntLaser(:,fsis(nanFind(i))));
        %determine if all nonlaser values are zero
        nolaserZeroFind = any(bigInt(:,fsis(nanFind(i))));
        if laserZeroFind == 0 && nolaserZeroFind == 0
            disp('No Spikes, Not adjusting')
        elseif laserZeroFind == 0 && nolaserZeroFind ~= 0
            disp('FULL SILENCING, ADJUSTING SLOPE')
            slopeStoreFSI(nanFind(i)) = 0;
            intStoreFSI(nanFind(i)) = 0;
        end
    end
end

%now we can generate an average based on the slopes and also based on all
%the points

%exclude negative slopes, as these are likely fucked. 
intStoreFSI(slopeStoreFSI<0) = NaN;
slopeStoreFSI(slopeStoreFSI < 0) = NaN;


avSlopeFSI = nanmean(slopeStoreFSI)
avIntFSI = nanmean(intStoreFSI)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);

hFig = figure;
hold on
for i = 1:length(fsis)
    plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.')
    plot(bigInt(:,fsis(i)),bigInt(:,fsis(i))*slopeStoreFSI(i) + intStoreFSI(i),'k')
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[avIntFSI maxVal*avSlopeFSI],'g')
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('FSI No Laser vs Laser')

spikeGraphName = 'peakBasedFSINoLaserVsLaserScatterWithAverageLine';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets look at average effects for MSNs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStore = [];
intStore = [];
valueStore=[];
for i = 1:length(msns)
    [b,bintr,bintjm] = gmregress(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,msns(i)),bigIntLaser(:,msns(i))];
    counter = counter + size(bigInt,1);
    slopeStore(i) = b(2);
    intStore(i) = b(1);
end

%now we can generate an average based on the slopes and also based on all
%the points

intStore(slopeStore<0) = NaN;
slopeStore(slopeStore < 0) = NaN;
intStore(slopeStore>3) = NaN;
slopeStore(slopeStore >3) = NaN;

avSlopeMSN = nanmean(slopeStore)
avIntMSN = nanmean(intStore)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);


hFig = figure;
hold on
for i = 1:length(msns)
    plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.')
    plot(bigInt(:,msns(i)),bigInt(:,msns(i))*slopeStore(i) + intStore(i),'k')
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[avIntMSN maxVal*avSlopeMSN],'g')
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('MSN No Laser vs Laser')


spikeGraphName = 'peakBasedMSNNoLaserVsLaserScatterWithAverageLine';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets plot out distribution of slopes. 
hFig = figure;
subplot(2,2,1)
hold on
hist(slopeStore,20)
tester = slopeStore;
tester(isnan(tester)) = [];
plot(nanmedian(tester),1,'g.')
plot(nanmean(tester),1.5,'r.')
title('MSN Slope Dist median green mean red')
subplot(2,2,2)
hold on
hist(intStore,20)
plot(nanmedian(intStore),1,'g.')
plot(nanmean(intStore),1.5,'r.')
title('MSN Int Dist')
subplot(2,2,3)
hold on
hist(slopeStoreFSI,20)
plot(nanmedian(slopeStoreFSI),1,'g.')
plot(nanmean(slopeStoreFSI),1.5,'r.')
title('FSI Slope Dist')
subplot(2,2,4)
hold on
hist(intStoreFSI,20)
plot(nanmedian(intStoreFSI),1,'g.')
plot(nanmean(intStoreFSI),1.5,'r.')
title('FSI Int Dist')

spikeGraphName = 'peakBasedMSNFSIDistributionSlopesAndIntercepts';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets try and do this thing where we find thresholds etc. First thing
%we need to do with this is clean up the data. 
fullSigTone(fullSigTone < 0.05) = 0;
fullSigTone(fullSigTone >0) = 1;
for i =1:length(msns)
    figure
    subplot(2,1,1)
    imagesc(fullSigTone(:,:,msns(i)))
    colormap('parula')
    title(unitName{i}, 'Interpreter', 'none')
    subplot(2,1,2)
    imagesc(fullSigLaser(:,:,msns(i)))
    colormap('parula')
end

%lets try and clean up this data, by eliminating things that arent
%connected to highest amplitude. This means going from the end, then
%marching down until I hit the first non-significant value. Remember that
%0s are significant, 1s are not

fullSigToneAnnot = fullSigTone;
fullSigLaserAnnot = fullSigLaser;

annotDBVals = zeros(length(fullSig),size(fullSigTone,1));
annotDBValsLaser = zeros(length(fullSig),size(fullSigTone,1));

for i = 1:length(fullSig)
    testData = fullSigTone(:,:,i);
    for j = 1:size(testData,1) %going by freq
        whileTrig = 0;
        whileCount = size(testData,2);
        while whileTrig == 0
            %determine if first value is zero or not
            if testData(j,whileCount) == 0
                disp(strcat('Significant Response at DB',num2str(whileCount)))
                if whileCount == 1
                    disp('Hit Lowest, Moving On to Next Freq')
                    annotDBVals(i,j) = whileCount;
                    whileTrig = 1;
                else
                    disp('Moving to Lower Amp')
                    whileCount = whileCount -1;
                end
            elseif testData(j,whileCount) == 1
                disp(strcat('First Non-Sig Resp at DB',num2str(whileCount)))
                annotDBVals(i,j) = whileCount+1;
                testData(j,1:whileCount) = 1;
                whileTrig = 1;
            end
        end
    end
    fullSigToneAnnot(:,:,i) = testData;
    disp('Finished a unit!')
end

for i = 1:length(fullSig)
    testData = fullSigLaser(:,:,i);
    for j = 1:size(testData,1) %going by freq
        whileTrig = 0;
        whileCount = size(testData,2);
        while whileTrig == 0
            %determine if first value is zero or not
            if testData(j,whileCount) == 0
                disp(strcat('Significant Response at DB',num2str(whileCount)))
                if whileCount == 1
                    disp('Hit Lowest, Moving On to Next Freq')
                    annotDBValsLaser(i,j) = whileCount;
                    whileTrig = 1;
                else
                    disp('Moving to Lower Amp')
                    whileCount = whileCount -1;
                end
            elseif testData(j,whileCount) == 1
                disp(strcat('First Non-Sig Resp at DB',num2str(whileCount)))
                annotDBValsLaser(i,j) = whileCount+1;
                testData(j,1:whileCount) = 1;
                whileTrig = 1;
            end
        end
    end
    fullSigLaserAnnot(:,:,i) = testData;
    disp('Finished a unit!')
end

for i =1:length(msns)
    figure
    subplot(2,1,1)
    imagesc(fullSigToneAnnot(:,:,msns(i)))
    colormap('parula')
    title(unitName{i}, 'Interpreter', 'none')
    subplot(2,1,2)
    imagesc(fullSigLaserAnnot(:,:,msns(i)))
    colormap('parula')
end

%this seems to work! Now lets eliminate units without enough responses
%anymore

for i = 1:length(fullSig)
    AnnotVals(i) = sum(sum(fullSigToneAnnot(:,:,i)));
end

figure
hist(AnnotVals,100)

%based on this, lets eliminate things with values greater than or equal to
%50

findBadAnnot = find(AnnotVals >=50);

annotMSNs = setdiff(msns,findBadAnnot);
annotFSIs = setdiff(fsis,findBadAnnot);

annotMSNDiffVals = annotDBVals(annotMSNs,:) -  annotDBValsLaser(annotMSNs,:);
annotFSIDiffVals = annotDBVals(annotFSIs,:) -  annotDBValsLaser(annotFSIs,:);

%calculate mean changes per unit
meanAnnotMSNDiff = mean(annotMSNDiffVals');
meanAnnotFSIDiff = mean(annotFSIDiffVals');

hFig = figure;
subplot(2,1,1)
hist(meanAnnotMSNDiff,10)
title('MSN Mean Difference in DB Threshold')
subplot(2,1,2)
hist(meanAnnotFSIDiff,10)
title('FSI Mean Difference In DB Change')

spikeGraphName = 'peakBasedFSIandMSNMeanDBThreshChange';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

signrank(meanAnnotMSNDiff)
signrank(meanAnnotFSIDiff)

clear


%% Now lets do it with the fast tone period?

testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
testNames = testNames(findString);

tonePeriod = 1;
numFiles = length(testNames);
counter = 1;
for j = 1:numFiles
    load(testNames{j})
    disp(strcat('Loading-',testNames{j}))
    numUnits = length(s.DesignationName);
    sigThresh = 0.05;
    valThresh = 8;
    sigValStore = [];
    for i = 1:numUnits
        %first, we need to see whether there is a significant response
        SigVals = s.(s.DesignationName{i}).BinSigVals(:,:,tonePeriod);
        %threshold it!
        SigVals(SigVals < sigThresh) = 0;
        SigVals(SigVals >0) = 1;
        SigValsLaser = s.(s.DesignationName{i}).BinSigValsLaser(:,:,tonePeriod);
        %threshold it!
        SigValsLaser(SigValsLaser < sigThresh) = 0;
        SigValsLaser(SigValsLaser >0) = 1;
        numSigVals = length(find(SigVals == 0));
        sigValStore(i) = numSigVals;
        %now lets use this as a cutoff
        if numSigVals > valThresh
            disp('Significant Crossing!')
            unitName{counter} = strcat(testNames{j},s.DesignationName{i});
            %now lets find the best frequency based on 70db response.
            [Y,I] = max(s.(s.DesignationName{i}).BinTone(:,9)); %I is best frequency
            if I
                intensityResponse(:,i) = s.(s.DesignationName{i}).BinDiff(I,:,tonePeriod);
                intensityResponseLaser(:,i) = s.(s.DesignationName{i}).BinDiffLaser(I,:,tonePeriod);
                bigInt(:,counter) = s.(s.DesignationName{i}).BinDiff(I,:,tonePeriod);
                bigIntLaser(:,counter) = s.(s.DesignationName{i}).BinDiffLaser(I,:,tonePeriod);
                bigSig(:,counter) = SigVals(I,:);
                bigSigLaser(:,counter) = SigValsLaser(I,:);
                fullSig(:,:,counter) = SigVals;
                fullSigTone(:,:,counter) = s.(s.DesignationName{i}).BinSigVals(:,:,2);
                fullSigLaser(:,:,counter) = SigValsLaser;
                fullVals(:,:,counter) = s.(s.DesignationName{i}).BinDiff(:,:,tonePeriod);
                fullValsLaser(:,:,counter) = s.(s.DesignationName{i}).BinDiffLaser(:,:,tonePeriod);
                cellType(counter) = masterData(i,7);
                counter = counter + 1;
                tempData = squeeze(s.(s.DesignationName{i}).BinSigVals(:,:,tonePeriod));
                minData = min(tempData);
                minData = find(minData < 0.05,1,'first');
                if minData
                    baseMin(counter) = minData;
                else
                    baseMin(counter) = 0;
                end
                tempData = squeeze(s.(s.DesignationName{i}).BinSigValsLaser(2:end,:,tonePeriod));
                minData = min(tempData);
                minData = find(minData < 0.05,1,'first');
                if minData
                    baseMinLaser(counter) = minData;
                else
                    baseMinLaser(counter) = 0;
                end
%                 figure
%                 subplot(2,1,1)
%                 plot(intensityResponse(:,i),'k')
%                 hold on
%                 plot(intensityResponseLaser(:,i) ,'g')
%                 
%                 subplot(2,1,2)
%                 plot(intensityResponseLaser(:,i) - intensityResponse(:,i))
            end
        end
    end

end

fsis = find(cellType == 1);
msns = find(cellType == 0);

hFig = figure;
subplot(2,1,1)
plot(bigInt(:,fsis))
title('FSI BF Amp Responses')
subplot(2,1,2)
plot(bigInt(:,msns))
title('MSN BF Amp Responses')

spikeGraphName = 'fastBinMSNFSIBFampResponses';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
subplot(2,1,1)
hold on
plot(mean(bigInt(:,fsis)'),'r')
plot(mean(bigInt(:,msns)'),'k')
plot(mean(bigIntLaser(:,msns)'),'g')
plot(mean(bigIntLaser(:,fsis)'),'m')
title('Mean Amp Resp MSNk FSIr')

subplot(2,1,2)
hold on
plot(mean(bigIntLaser(:,fsis)')-mean(bigInt(:,fsis)'),'r')
plot(mean(bigIntLaser(:,msns)')-mean(bigInt(:,msns)'),'k')
title('Mean Diff Amp Resp MSNk FSIr')

spikeGraphName = 'fastBinMSNFSIPopAvAmpResponses';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
subplot(2,1,1)
hold on
for i = 1:length(fsis)
    plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('Ind FSI Amp Resp Scatter Laser vs NoLaser')
xlabel('spikes no laser')
ylabel('spikes with laser')

figure
hold on
for i = 1:length(msns)
    plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('Ind MSN Amp Resp Scatter Laser vs NoLaser')
xlabel('spikes no laser')
ylabel('spikes with laser')


spikeGraphName = 'fastBinMSNFSIIndivSpikeScatterLaserNoLaser';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
%lets look at average effects for FSIs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStoreFSI = [];
intStoreFSI = [];
valueStore=[];
for i = 1:length(fsis)
    [b,bintr,bintjm] = gmregress(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,fsis(i)),bigIntLaser(:,fsis(i))];
    counter = counter + size(bigInt,1);
    slopeStoreFSI(i) = b(2);
    intStoreFSI(i) = b(1);
end

%determine if there are any values where slope is NaN
nanFind = find(isnan(slopeStoreFSI));
if nanFind
    disp('Some GMREGRESS FAIL FSIs')
    for i = 1:length(nanFind)
        %determine if all laser values are zero
        laserZeroFind = any(bigIntLaser(:,fsis(nanFind(i))));
        %determine if all nonlaser values are zero
        nolaserZeroFind = any(bigInt(:,fsis(nanFind(i))));
        if laserZeroFind == 0 && nolaserZeroFind == 0
            disp('No Spikes, Not adjusting')
        elseif laserZeroFind == 0 && nolaserZeroFind ~= 0
            disp('FULL SILENCING, ADJUSTING SLOPE')
            slopeStoreFSI(nanFind(i)) = 0;
            intStoreFSI(nanFind(i)) = 0;
        end
    end
end

%now we can generate an average based on the slopes and also based on all
%the points

%exclude negative slopes, as these are likely fucked. 
intStoreFSI(slopeStoreFSI<0) = NaN;
slopeStoreFSI(slopeStoreFSI < 0) = NaN;


avSlopeFSI = nanmean(slopeStoreFSI)
avIntFSI = nanmean(intStoreFSI)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);

hFig = figure;
hold on
for i = 1:length(fsis)
    plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.')
    plot(bigInt(:,fsis(i)),bigInt(:,fsis(i))*slopeStoreFSI(i) + intStoreFSI(i),'k')
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[avIntFSI maxVal*avSlopeFSI],'g')
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('FSI No Laser vs Laser')

spikeGraphName = 'fastBinFSINoLaserVsLaserScatterWithAverageLine';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets look at average effects for MSNs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStore = [];
intStore = [];
valueStore=[];
for i = 1:length(msns)
    [b,bintr,bintjm] = gmregress(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,msns(i)),bigIntLaser(:,msns(i))];
    counter = counter + size(bigInt,1);
    slopeStore(i) = b(2);
    intStore(i) = b(1);
end

%now we can generate an average based on the slopes and also based on all
%the points

intStore(slopeStore<0) = NaN;
slopeStore(slopeStore < 0) = NaN;
intStore(slopeStore>3) = NaN;
slopeStore(slopeStore >3) = NaN;

avSlopeMSN = nanmean(slopeStore)
avIntMSN = nanmean(intStore)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);


hFig = figure;
hold on
for i = 1:length(msns)
    plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.')
    plot(bigInt(:,msns(i)),bigInt(:,msns(i))*slopeStore(i) + intStore(i),'k')
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[avIntMSN maxVal*avSlopeMSN],'g')
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('MSN No Laser vs Laser')


spikeGraphName = 'fastBinMSNNoLaserVsLaserScatterWithAverageLine';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets plot out distribution of slopes. 
hFig = figure;
subplot(2,2,1)
hold on
hist(slopeStore,20)
tester = slopeStore;
tester(isnan(tester)) = [];
plot(nanmedian(tester),1,'g.')
plot(nanmean(tester),1.5,'r.')
title('MSN Slope Dist median green mean red')
subplot(2,2,2)
hold on
hist(intStore,20)
plot(nanmedian(intStore),1,'g.')
plot(nanmean(intStore),1.5,'r.')
title('MSN Int Dist')
subplot(2,2,3)
hold on
hist(slopeStoreFSI,20)
plot(nanmedian(slopeStoreFSI),1,'g.')
plot(nanmean(slopeStoreFSI),1.5,'r.')
title('FSI Slope Dist')
subplot(2,2,4)
hold on
hist(intStoreFSI,20)
plot(nanmedian(intStoreFSI),1,'g.')
plot(nanmean(intStoreFSI),1.5,'r.')
title('FSI Int Dist')

spikeGraphName = 'fastBinMSNFSIDistributionSlopesAndIntercepts';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets try and do this thing where we find thresholds etc. First thing
%we need to do with this is clean up the data. 
fullSigTone(fullSigTone < 0.05) = 0;
fullSigTone(fullSigTone >0) = 1;
for i =1:length(msns)
    figure
    subplot(2,1,1)
    imagesc(fullSigTone(:,:,msns(i)))
    colormap('parula')
    title(unitName{i}, 'Interpreter', 'none')
    subplot(2,1,2)
    imagesc(fullSigLaser(:,:,msns(i)))
    colormap('parula')
end

%lets try and clean up this data, by eliminating things that arent
%connected to highest amplitude. This means going from the end, then
%marching down until I hit the first non-significant value. Remember that
%0s are significant, 1s are not

fullSigToneAnnot = fullSigTone;
fullSigLaserAnnot = fullSigLaser;

annotDBVals = zeros(length(fullSig),size(fullSigTone,1));
annotDBValsLaser = zeros(length(fullSig),size(fullSigTone,1));

for i = 1:length(fullSig)
    testData = fullSigTone(:,:,i);
    for j = 1:size(testData,1) %going by freq
        whileTrig = 0;
        whileCount = size(testData,2);
        while whileTrig == 0
            %determine if first value is zero or not
            if testData(j,whileCount) == 0
                disp(strcat('Significant Response at DB',num2str(whileCount)))
                if whileCount == 1
                    disp('Hit Lowest, Moving On to Next Freq')
                    annotDBVals(i,j) = whileCount;
                    whileTrig = 1;
                else
                    disp('Moving to Lower Amp')
                    whileCount = whileCount -1;
                end
            elseif testData(j,whileCount) == 1
                disp(strcat('First Non-Sig Resp at DB',num2str(whileCount)))
                annotDBVals(i,j) = whileCount+1;
                testData(j,1:whileCount) = 1;
                whileTrig = 1;
            end
        end
    end
    fullSigToneAnnot(:,:,i) = testData;
    disp('Finished a unit!')
end

for i = 1:length(fullSig)
    testData = fullSigLaser(:,:,i);
    for j = 1:size(testData,1) %going by freq
        whileTrig = 0;
        whileCount = size(testData,2);
        while whileTrig == 0
            %determine if first value is zero or not
            if testData(j,whileCount) == 0
                disp(strcat('Significant Response at DB',num2str(whileCount)))
                if whileCount == 1
                    disp('Hit Lowest, Moving On to Next Freq')
                    annotDBValsLaser(i,j) = whileCount;
                    whileTrig = 1;
                else
                    disp('Moving to Lower Amp')
                    whileCount = whileCount -1;
                end
            elseif testData(j,whileCount) == 1
                disp(strcat('First Non-Sig Resp at DB',num2str(whileCount)))
                annotDBValsLaser(i,j) = whileCount+1;
                testData(j,1:whileCount) = 1;
                whileTrig = 1;
            end
        end
    end
    fullSigLaserAnnot(:,:,i) = testData;
    disp('Finished a unit!')
end

for i =1:length(msns)
    figure
    subplot(2,1,1)
    imagesc(fullSigToneAnnot(:,:,msns(i)))
    colormap('parula')
    title(unitName{i}, 'Interpreter', 'none')
    subplot(2,1,2)
    imagesc(fullSigLaserAnnot(:,:,msns(i)))
    colormap('parula')
end

%this seems to work! Now lets eliminate units without enough responses
%anymore

for i = 1:length(fullSig)
    AnnotVals(i) = sum(sum(fullSigToneAnnot(:,:,i)));
end

figure
hist(AnnotVals,100)

%based on this, lets eliminate things with values greater than or equal to
%50

findBadAnnot = find(AnnotVals >=50);

annotMSNs = setdiff(msns,findBadAnnot);
annotFSIs = setdiff(fsis,findBadAnnot);

annotMSNDiffVals = annotDBVals(annotMSNs,:) -  annotDBValsLaser(annotMSNs,:);
annotFSIDiffVals = annotDBVals(annotFSIs,:) -  annotDBValsLaser(annotFSIs,:);

%calculate mean changes per unit
meanAnnotMSNDiff = mean(annotMSNDiffVals');
meanAnnotFSIDiff = mean(annotFSIDiffVals');

hFig = figure;
subplot(2,1,1)
hist(meanAnnotMSNDiff,10)
title('MSN Mean Difference in DB Threshold')
subplot(2,1,2)
hist(meanAnnotFSIDiff,10)
title('FSI Mean Difference In DB Change')

spikeGraphName = 'fastBinFSIandMSNMeanDBThreshChange';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

signrank(meanAnnotMSNDiff)
signrank(meanAnnotFSIDiff)

clear

%% Now plot for general period
testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
testNames = testNames(findString);

tonePeriod = 3;
numFiles = length(testNames);
counter = 1;
for j = 1:numFiles
    load(testNames{j})
    disp(strcat('Loading-',testNames{j}))
    numUnits = length(s.DesignationName);
    sigThresh = 0.05;
    valThresh = 8;
    sigValStore = [];
    for i = 1:numUnits
        %first, we need to see whether there is a significant response
        SigVals = s.(s.DesignationName{i}).BinSigVals(:,:,tonePeriod);
        %threshold it!
        SigVals(SigVals < sigThresh) = 0;
        SigVals(SigVals >0) = 1;
        SigValsLaser = s.(s.DesignationName{i}).BinSigValsLaser(:,:,tonePeriod);
        %threshold it!
        SigValsLaser(SigValsLaser < sigThresh) = 0;
        SigValsLaser(SigValsLaser >0) = 1;
        numSigVals = length(find(SigVals == 0));
        sigValStore(i) = numSigVals;
        %now lets use this as a cutoff
        if numSigVals > valThresh
            disp('Significant Crossing!')
            unitName{counter} = strcat(testNames{j},s.DesignationName{i});
            %now lets find the best frequency based on 70db response.
            [Y,I] = max(s.(s.DesignationName{i}).BinTone(:,9)); %I is best frequency
            if I
                intensityResponse(:,i) = s.(s.DesignationName{i}).BinDiff(I,:,tonePeriod);
                intensityResponseLaser(:,i) = s.(s.DesignationName{i}).BinDiffLaser(I,:,tonePeriod);
                bigInt(:,counter) = s.(s.DesignationName{i}).BinDiff(I,:,tonePeriod);
                bigIntLaser(:,counter) = s.(s.DesignationName{i}).BinDiffLaser(I,:,tonePeriod);
                bigSig(:,counter) = SigVals(I,:);
                bigSigLaser(:,counter) = SigValsLaser(I,:);
                fullSig(:,:,counter) = SigVals;
                fullSigTone(:,:,counter) = s.(s.DesignationName{i}).BinSigVals(:,:,2);
                fullSigLaser(:,:,counter) = SigValsLaser;
                fullVals(:,:,counter) = s.(s.DesignationName{i}).BinDiff(:,:,tonePeriod);
                fullValsLaser(:,:,counter) = s.(s.DesignationName{i}).BinDiffLaser(:,:,tonePeriod);
                cellType(counter) = masterData(i,7);
                counter = counter + 1;
                tempData = squeeze(s.(s.DesignationName{i}).BinSigVals(:,:,tonePeriod));
                minData = min(tempData);
                minData = find(minData < 0.05,1,'first');
                if minData
                    baseMin(counter) = minData;
                else
                    baseMin(counter) = 0;
                end
                tempData = squeeze(s.(s.DesignationName{i}).BinSigValsLaser(2:end,:,tonePeriod));
                minData = min(tempData);
                minData = find(minData < 0.05,1,'first');
                if minData
                    baseMinLaser(counter) = minData;
                else
                    baseMinLaser(counter) = 0;
                end
%                 figure
%                 subplot(2,1,1)
%                 plot(intensityResponse(:,i),'k')
%                 hold on
%                 plot(intensityResponseLaser(:,i) ,'g')
%                 
%                 subplot(2,1,2)
%                 plot(intensityResponseLaser(:,i) - intensityResponse(:,i))
            end
        end
    end

end

fsis = find(cellType == 1);
msns = find(cellType == 0);

hFig = figure;
subplot(2,1,1)
plot(bigInt(:,fsis))
title('FSI BF Amp Responses')
subplot(2,1,2)
plot(bigInt(:,msns))
title('MSN BF Amp Responses')

spikeGraphName = 'genBinMSNFSIBFampResponses';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
subplot(2,1,1)
hold on
plot(mean(bigInt(:,fsis)'),'r')
plot(mean(bigInt(:,msns)'),'k')
plot(mean(bigIntLaser(:,msns)'),'g')
plot(mean(bigIntLaser(:,fsis)'),'m')
title('Mean Amp Resp MSNk FSIr')

subplot(2,1,2)
hold on
plot(mean(bigIntLaser(:,fsis)')-mean(bigInt(:,fsis)'),'r')
plot(mean(bigIntLaser(:,msns)')-mean(bigInt(:,msns)'),'k')
title('Mean Diff Amp Resp MSNk FSIr')

spikeGraphName = 'genBinMSNFSIPopAvAmpResponses';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

hFig = figure;
subplot(2,1,1)
hold on
for i = 1:length(fsis)
    plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('Ind FSI Amp Resp Scatter Laser vs NoLaser')
xlabel('spikes no laser')
ylabel('spikes with laser')

figure
hold on
for i = 1:length(msns)
    plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('Ind MSN Amp Resp Scatter Laser vs NoLaser')
xlabel('spikes no laser')
ylabel('spikes with laser')


spikeGraphName = 'genBinMSNFSIIndivSpikeScatterLaserNoLaser';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')
%lets look at average effects for FSIs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStoreFSI = [];
intStoreFSI = [];
valueStore=[];
for i = 1:length(fsis)
    [b,bintr,bintjm] = gmregress(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,fsis(i)),bigIntLaser(:,fsis(i))];
    counter = counter + size(bigInt,1);
    slopeStoreFSI(i) = b(2);
    intStoreFSI(i) = b(1);
end

%determine if there are any values where slope is NaN
nanFind = find(isnan(slopeStoreFSI));
if nanFind
    disp('Some GMREGRESS FAIL FSIs')
    for i = 1:length(nanFind)
        %determine if all laser values are zero
        laserZeroFind = any(bigIntLaser(:,fsis(nanFind(i))));
        %determine if all nonlaser values are zero
        nolaserZeroFind = any(bigInt(:,fsis(nanFind(i))));
        if laserZeroFind == 0 && nolaserZeroFind == 0
            disp('No Spikes, Not adjusting')
        elseif laserZeroFind == 0 && nolaserZeroFind ~= 0
            disp('FULL SILENCING, ADJUSTING SLOPE')
            slopeStoreFSI(nanFind(i)) = 0;
            intStoreFSI(nanFind(i)) = 0;
        end
    end
end

%now we can generate an average based on the slopes and also based on all
%the points

%exclude negative slopes, as these are likely fucked. 
intStoreFSI(slopeStoreFSI<0) = NaN;
slopeStoreFSI(slopeStoreFSI < 0) = NaN;


avSlopeFSI = nanmean(slopeStoreFSI)
avIntFSI = nanmean(intStoreFSI)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);

hFig = figure;
hold on
for i = 1:length(fsis)
    plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.')
    plot(bigInt(:,fsis(i)),bigInt(:,fsis(i))*slopeStoreFSI(i) + intStoreFSI(i),'k')
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[avIntFSI maxVal*avSlopeFSI],'g')
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('FSI No Laser vs Laser')

spikeGraphName = 'genBinFSINoLaserVsLaserScatterWithAverageLine';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%lets look at average effects for MSNs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStore = [];
intStore = [];
valueStore=[];
for i = 1:length(msns)
    [b,bintr,bintjm] = gmregress(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,msns(i)),bigIntLaser(:,msns(i))];
    counter = counter + size(bigInt,1);
    slopeStore(i) = b(2);
    intStore(i) = b(1);
end

%now we can generate an average based on the slopes and also based on all
%the points

intStore(slopeStore<0) = NaN;
slopeStore(slopeStore < 0) = NaN;
intStore(slopeStore>3) = NaN;
slopeStore(slopeStore >3) = NaN;

avSlopeMSN = nanmean(slopeStore)
avIntMSN = nanmean(intStore)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);


hFig = figure;
hold on
for i = 1:length(msns)
    plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.')
    plot(bigInt(:,msns(i)),bigInt(:,msns(i))*slopeStore(i) + intStore(i),'k')
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[avIntMSN maxVal*avSlopeMSN],'g')
plot([0 maxVal],[0 maxVal],'r')
axis equal
xlim([0 maxVal])
ylim([0 maxVal])
title('MSN No Laser vs Laser')


spikeGraphName = 'genBinMSNNoLaserVsLaserScatterWithAverageLine';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets plot out distribution of slopes. 
hFig = figure;
subplot(2,2,1)
hold on
hist(slopeStore,20)
tester = slopeStore;
tester(isnan(tester)) = [];
plot(nanmedian(tester),1,'g.')
plot(nanmean(tester),1.5,'r.')
title('MSN Slope Dist median green mean red')
subplot(2,2,2)
hold on
hist(intStore,20)
plot(nanmedian(intStore),1,'g.')
plot(nanmean(intStore),1.5,'r.')
title('MSN Int Dist')
subplot(2,2,3)
hold on
hist(slopeStoreFSI,20)
plot(nanmedian(slopeStoreFSI),1,'g.')
plot(nanmean(slopeStoreFSI),1.5,'r.')
title('FSI Slope Dist')
subplot(2,2,4)
hold on
hist(intStoreFSI,20)
plot(nanmedian(intStoreFSI),1,'g.')
plot(nanmean(intStoreFSI),1.5,'r.')
title('FSI Int Dist')

spikeGraphName = 'genBinMSNFSIDistributionSlopesAndIntercepts';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

%now lets try and do this thing where we find thresholds etc. First thing
%we need to do with this is clean up the data. 
fullSigTone(fullSigTone < 0.05) = 0;
fullSigTone(fullSigTone >0) = 1;
for i =1:length(msns)
    figure
    subplot(2,1,1)
    imagesc(fullSigTone(:,:,msns(i)))
    colormap('parula')
    title(unitName{i}, 'Interpreter', 'none')
    subplot(2,1,2)
    imagesc(fullSigLaser(:,:,msns(i)))
    colormap('parula')
end

%lets try and clean up this data, by eliminating things that arent
%connected to highest amplitude. This means going from the end, then
%marching down until I hit the first non-significant value. Remember that
%0s are significant, 1s are not

fullSigToneAnnot = fullSigTone;
fullSigLaserAnnot = fullSigLaser;

annotDBVals = zeros(length(fullSig),size(fullSigTone,1));
annotDBValsLaser = zeros(length(fullSig),size(fullSigTone,1));

for i = 1:length(fullSig)
    testData = fullSigTone(:,:,i);
    for j = 1:size(testData,1) %going by freq
        whileTrig = 0;
        whileCount = size(testData,2);
        while whileTrig == 0
            %determine if first value is zero or not
            if testData(j,whileCount) == 0
                disp(strcat('Significant Response at DB',num2str(whileCount)))
                if whileCount == 1
                    disp('Hit Lowest, Moving On to Next Freq')
                    annotDBVals(i,j) = whileCount;
                    whileTrig = 1;
                else
                    disp('Moving to Lower Amp')
                    whileCount = whileCount -1;
                end
            elseif testData(j,whileCount) == 1
                disp(strcat('First Non-Sig Resp at DB',num2str(whileCount)))
                annotDBVals(i,j) = whileCount+1;
                testData(j,1:whileCount) = 1;
                whileTrig = 1;
            end
        end
    end
    fullSigToneAnnot(:,:,i) = testData;
    disp('Finished a unit!')
end

for i = 1:length(fullSig)
    testData = fullSigLaser(:,:,i);
    for j = 1:size(testData,1) %going by freq
        whileTrig = 0;
        whileCount = size(testData,2);
        while whileTrig == 0
            %determine if first value is zero or not
            if testData(j,whileCount) == 0
                disp(strcat('Significant Response at DB',num2str(whileCount)))
                if whileCount == 1
                    disp('Hit Lowest, Moving On to Next Freq')
                    annotDBValsLaser(i,j) = whileCount;
                    whileTrig = 1;
                else
                    disp('Moving to Lower Amp')
                    whileCount = whileCount -1;
                end
            elseif testData(j,whileCount) == 1
                disp(strcat('First Non-Sig Resp at DB',num2str(whileCount)))
                annotDBValsLaser(i,j) = whileCount+1;
                testData(j,1:whileCount) = 1;
                whileTrig = 1;
            end
        end
    end
    fullSigLaserAnnot(:,:,i) = testData;
    disp('Finished a unit!')
end

for i =1:length(msns)
    figure
    subplot(2,1,1)
    imagesc(fullSigToneAnnot(:,:,msns(i)))
    colormap('parula')
    title(unitName{i}, 'Interpreter', 'none')
    subplot(2,1,2)
    imagesc(fullSigLaserAnnot(:,:,msns(i)))
    colormap('parula')
end

%this seems to work! Now lets eliminate units without enough responses
%anymore

for i = 1:length(fullSig)
    AnnotVals(i) = sum(sum(fullSigToneAnnot(:,:,i)));
end

figure
hist(AnnotVals,100)

%based on this, lets eliminate things with values greater than or equal to
%50

findBadAnnot = find(AnnotVals >=50);

annotMSNs = setdiff(msns,findBadAnnot);
annotFSIs = setdiff(fsis,findBadAnnot);

annotMSNDiffVals = annotDBVals(annotMSNs,:) -  annotDBValsLaser(annotMSNs,:);
annotFSIDiffVals = annotDBVals(annotFSIs,:) -  annotDBValsLaser(annotFSIs,:);

%calculate mean changes per unit
meanAnnotMSNDiff = mean(annotMSNDiffVals');
meanAnnotFSIDiff = mean(annotFSIDiffVals');

hFig = figure;
subplot(2,1,1)
hist(meanAnnotMSNDiff,10)
title('MSN Mean Difference in DB Threshold')
subplot(2,1,2)
hist(meanAnnotFSIDiff,10)
title('FSI Mean Difference In DB Change')

spikeGraphName = 'genBinFSIandMSNMeanDBThreshChange';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')

signrank(meanAnnotMSNDiff)
signrank(meanAnnotFSIDiff)

clear


%% now lets look at some individual files. 
load('181114_ML180820I_L_AudStr_pen2_3320_3mWPVHaloTuningAltLaserFullTuningAnalysis.mat')

%lets start with 13 cluster 2
baseline = s.nt13cluster2.BinDiff(:,:,2);
laser = s.nt13cluster2.BinDiffLaser(:,:,2);
baselineSig = s.nt13cluster2.BinSigVals(:,:,2);
laserSig = s.nt13cluster2.BinSigValsLaser(:,:,2);
%lets threshold to some value
tarSig = 0.05;
baselineSig(baselineSig < tarSig) = 0;
baselineSig(baselineSig > 0) = 1;
laserSig(laserSig < tarSig) = 0;
laserSig(laserSig > 0) = 1;
maxval = max([max(max(baseline)),max(max(laser))]);
minval = min([min(min(baseline)),min(min(laser))]);

figure
%plot baseline subtracted values
subplot(4,2,1)
imagesc(baseline,[minval,maxval])
colormap('parula')
colorbar
title('BaseSub NoLaser')
subplot(4,2,3)
imagesc(laser,[minval,maxval])
colormap('parula')
colorbar
title('BaseSub WithLaser')

for i = 1:6
    subplot(12,2,11+2*i)
    hold on
    plot(baseline(i,:),'k.-')
    %now find significant responses
    testFind = find(baselineSig(i,:) == 0);
    plot(testFind,baseline(i,testFind),'k*')
    plot(laser(i,:),'g.-')
    %now find significant responses
    testFind = find(laserSig(i,:) == 0);
    plot(testFind,laser(i,testFind),'g*')
end
%now plot non-baseline subtracted
baseline = s.nt13cluster2.BinTone(:,:);
laser = s.nt13cluster2.BinToneLaser(:,:);
baselineSig = s.nt13cluster2.BinSigVals(:,:,2);
laserSig = s.nt13cluster2.BinSigValsLaser(:,:,2);
%lets threshold to some value
tarSig = 0.05;
baselineSig(baselineSig < tarSig) = 0;
baselineSig(baselineSig > 0) = 1;
laserSig(laserSig < tarSig) = 0;
laserSig(laserSig > 0) = 1;
maxval = max([max(max(baseline)),max(max(laser))]);
minval = min([min(min(baseline)),min(min(laser))]);

subplot(4,2,2)
imagesc(baseline,[minval,maxval])
colormap('parula')
colorbar
title('NoBaseSub NoLaser')
subplot(4,2,4)
imagesc(laser,[minval,maxval])
colormap('parula')
colorbar
title('NoBaseSub WithLaser')

for i = 1:6
    subplot(12,2,12+2*i)
    hold on
    plot(baseline(i,:),'k.-')
    %now find significant responses
    testFind = find(baselineSig(i,:) == 0);
    plot(testFind,baseline(i,testFind),'k*')
    plot(laser(i,:),'g.-')
    %now find significant responses
    testFind = find(laserSig(i,:) == 0);
    plot(testFind,laser(i,testFind),'g*')
end

%now 14.1
baseline = s.nt14cluster1.BinDiff(:,:,2);
laser = s.nt14cluster1.BinDiffLaser(:,:,2);
baselineSig = s.nt14cluster1.BinSigVals(:,:,2);
laserSig = s.nt14cluster1.BinSigValsLaser(:,:,2);
%lets threshold to some value
tarSig = 0.05;
baselineSig(baselineSig < tarSig) = 0;
baselineSig(baselineSig > 0) = 1;
laserSig(laserSig < tarSig) = 0;
laserSig(laserSig > 0) = 1;
maxval = max([max(max(baseline)),max(max(laser))]);
minval = min([min(min(baseline)),min(min(laser))]);

figure
subplot(2,1,1)
imagesc(baseline,[minval,maxval])
colormap('parula')
colorbar
subplot(2,1,2)
imagesc(laser,[minval,maxval])
colormap('parula')
colorbar

figure
for i = 1:6
    subplot(6,1,i)
    hold on
    plot(baseline(i,:),'k.-')
    %now find significant responses
    testFind = find(baselineSig(i,:) == 0);
    plot(testFind,baseline(i,testFind),'k*')
    plot(laser(i,:),'g.-')
    %now find significant responses
    testFind = find(laserSig(i,:) == 0);
    plot(testFind,laser(i,testFind),'g*')
end


%now 15.1
baseline = s.nt15cluster1.BinDiff(:,:,2);
laser = s.nt15cluster1.BinDiffLaser(:,:,2);
baselineSig = s.nt15cluster1.BinSigVals(:,:,2);
laserSig = s.nt15cluster1.BinSigValsLaser(:,:,2);
%lets threshold to some value
tarSig = 0.05;
baselineSig(baselineSig < tarSig) = 0;
baselineSig(baselineSig > 0) = 1;
laserSig(laserSig < tarSig) = 0;
laserSig(laserSig > 0) = 1;
maxval = max([max(max(baseline)),max(max(laser))]);
minval = min([min(min(baseline)),min(min(laser))]);

figure
subplot(2,1,1)
imagesc(baseline,[minval,maxval])
colormap('parula')
colorbar
subplot(2,1,2)
imagesc(laser,[minval,maxval])
colormap('parula')
colorbar

figure
for i = 1:6
    subplot(6,1,i)
    hold on
    plot(baseline(i,:),'k.-')
    %now find significant responses
    testFind = find(baselineSig(i,:) == 0);
    plot(testFind,baseline(i,testFind),'k*')
    plot(laser(i,:),'g.-')
    %now find significant responses
    testFind = find(laserSig(i,:) == 0);
    plot(testFind,laser(i,testFind),'g*')
end


%now 16.4
baseline = s.nt16cluster4.BinDiff(:,:,2);
laser = s.nt16cluster4.BinDiffLaser(:,:,2);
baselineSig = s.nt16cluster4.BinSigVals(:,:,2);
laserSig = s.nt16cluster4.BinSigValsLaser(:,:,2);
%lets threshold to some value
tarSig = 0.05;
baselineSig(baselineSig < tarSig) = 0;
baselineSig(baselineSig > 0) = 1;
laserSig(laserSig < tarSig) = 0;
laserSig(laserSig > 0) = 1;
maxval = max([max(max(baseline)),max(max(laser))]);
minval = min([min(min(baseline)),min(min(laser))]);

figure
subplot(2,1,1)
imagesc(baseline,[minval,maxval])
colormap('parula')
colorbar
subplot(2,1,2)
imagesc(laser,[minval,maxval])
colormap('parula')
colorbar

figure
for i = 1:6
    subplot(6,1,i)
    hold on
    plot(baseline(i,:),'k.-')
    %now find significant responses
    testFind = find(baselineSig(i,:) == 0);
    plot(testFind,baseline(i,testFind),'k*')
    plot(laser(i,:),'g.-')
    %now find significant responses
    testFind = find(laserSig(i,:) == 0);
    plot(testFind,laser(i,testFind),'g*')
end


%based on these data, it seems like i'll have a hard time making a case
%that there is a creep of the amplitude threshold. See this on rare
%occasion, but doesnt seem to be a common pattern. 








%% Lets make a version to apply to baseline recording data.

testNames = what;
testNames = testNames.mat;
[findString] = functionCellStringFind(testNames,'Analysis');
testNames = testNames(findString);


numFiles = length(testNames);
counter = 1;
for j = 1:numFiles
    load(testNames{j})
    disp(strcat('Loading-',testNames{j}))
    numUnits = length(s.DesignationName);
    sigThresh = 0.05;
    valThresh = 8;
    sigValStore = [];
    for i = 1:numUnits
        %first, we need to see whether there is a significant response
        SigVals = s.(s.DesignationName{i}).BinSigVals(:,:,2);
        %threshold it!
        SigVals(SigVals < sigThresh) = 0;
        SigVals(SigVals >0) = 1;
        numSigVals = length(find(SigVals == 0));
        sigValStore(i) = numSigVals;
        %now lets use this as a cutoff
        if numSigVals > valThresh
            disp('Significant Crossing!')
            %now lets find the best frequency based on 70db response.
            [Y,I] = max(s.(s.DesignationName{i}).BinTone(:,end)); %I is best frequency
            if I
                intensityResponse(:,i) = s.(s.DesignationName{i}).BinTone(I,:);
%                 intensityResponseLaser(:,i) = s.(s.DesignationName{i}).BinToneLaser(I,:);
                bigInt(:,counter) = s.(s.DesignationName{i}).BinTone(I,:);
%                 bigIntLaser(:,counter) = s.(s.DesignationName{i}).BinToneLaser(I,:);
                cellType(counter) = masterData(i,7);
                counter = counter + 1;
%                 figure
%                 subplot(2,1,1)
%                 plot(intensityResponse(:,i),'k')
%                 hold on
%                 plot(intensityResponseLaser(:,i) ,'g')
%                 
%                 subplot(2,1,2)
%                 plot(intensityResponseLaser(:,i) - intensityResponse(:,i))
            end
        end
    end

end


fsis = find(cellType == 1);
msns = find(cellType == 0);


figure
subplot(2,1,1)
plot(bigInt(:,fsis))
title('FSI BF Amp Responses')
subplot(2,1,2)
plot(bigInt(:,msns))
title('MSN BF Amp Responses')

figure
hold on
plot(mean(bigInt(:,fsis)'),'r')
plot(mean(bigInt(:,msns)'),'k')
title('Mean Amp Resp MSNk FSIr')
