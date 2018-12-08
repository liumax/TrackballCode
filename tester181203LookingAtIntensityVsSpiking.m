%lets look at overall plots! 
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
        SigValsLaser = s.(s.DesignationName{i}).BinSigVals(:,:,tonePeriod);
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
plot(mean(bigIntLaser(:,msns)'),'g')
plot(mean(bigIntLaser(:,fsis)'),'m')
title('Mean Amp Resp MSNk FSIr')

figure
hold on
plot(mean(bigIntLaser(:,fsis)')-mean(bigInt(:,fsis)'),'r')
plot(mean(bigIntLaser(:,msns)')-mean(bigInt(:,msns)'),'k')
title('Mean Diff Amp Resp MSNk FSIr')

figure
hold on
for i = 1:length(fsis)
    plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[0 maxVal],'r')

figure
hold on
for i = 1:length(msns)
    plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[0 maxVal],'r')

%lets look at average effects for FSIs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStore = [];
intStore = [];
valueStore=[];
for i = 1:length(fsis)
    [b,bintr,bintjm] = gmregress(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,fsis(i)),bigIntLaser(:,fsis(i))];
    counter = counter + size(bigInt,1);
    slopeStore(i) = b(2);
    intStore(i) = b(1);
end

%now we can generate an average based on the slopes and also based on all
%the points

avSlopeFSI = nanmean(slopeStore)
avIntFSI = nanmean(intStore)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);

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

avSlopeMSN = nanmean(slopeStore)
avIntMSN = nanmean(intStore)
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);


%now lets try and do this thing where we find thresholds etc. First thing
%we need to do with this is clean up the data. 
fullSigTone(fullSigTone < 0.05) = 0;
fullSigTone(fullSigTone >0) = 1;
for i =1:length(fullSig)
    figure
    subplot(2,1,1)
    imagesc(fullSig(:,:,i))
    colormap('parula')
    title(unitName{i}, 'Interpreter', 'none')
    subplot(2,1,2)
    imagesc(fullSigTone(:,:,i))
    colormap('parula')
end



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
plot(mean(bigIntLaser(:,msns)'),'g')
plot(mean(bigIntLaser(:,fsis)'),'m')
title('Mean Amp Resp MSNk FSIr')

figure
hold on
plot(mean(bigIntLaser(:,fsis)')-mean(bigInt(:,fsis)'),'r')
plot(mean(bigIntLaser(:,msns)')-mean(bigInt(:,msns)'),'k')
title('Mean Diff Amp Resp MSNk FSIr')

figure
hold on
for i = 1:length(fsis)
    plot(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,fsis)));
plot([0 maxVal],[0 maxVal],'r')

figure
hold on
for i = 1:length(msns)
    plot(bigInt(:,msns(i)),bigIntLaser(:,msns(i)),'k.-')
end
maxVal = max(max(bigIntLaser(:,msns)));
plot([0 maxVal],[0 maxVal],'r')

%lets look at average effects for FSIs

% [b,bintr,bintjm] = gmregress(valStore(:,1),valStore(:,2),sigVal);
sigVal = 0.05;
counter = 1;
slopeStore = [];
intStore = [];
valueStore=[];
for i = 1:length(fsis)
    [b,bintr,bintjm] = gmregress(bigInt(:,fsis(i)),bigIntLaser(:,fsis(i)),sigVal);
    %lets also store the values.
    valueStore(counter:counter + size(bigInt,1) - 1,:) = [bigInt(:,fsis(i)),bigIntLaser(:,fsis(i))];
    counter = counter + size(bigInt,1);
    slopeStore(i) = b(2);
    intStore(i) = b(1);
end

%now we can generate an average based on the slopes and also based on all
%the points

avSlope = nanmean(slopeStore);
avInt = nanmean(intStore);
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);

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

avSlope = nanmean(slopeStore);
avInt = nanmean(intStore);
[b,bintr,bintjm] = gmregress(valueStore(:,1),valueStore(:,2),sigVal);



%% now lets look at some individual files. 
load('181114_ML180820I_L_AudStr_pen2_3320_3mWPVHaloTuningAltLaserFullTuningAnalysis (1).mat')

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
