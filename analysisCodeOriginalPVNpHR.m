%first navigate to correct folder
cd /Users/maxliu/181212LookingAtOriginalPVNPHR

%now lets go through the files and extract data. This should be updated
%data with the correct system for finding significant responses

targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);


numFiles = length(targetFiles);

%okay, now that we have files, lets do something with them! 


%first thing i'm interested in checking is whether I can make significant
%regressions using the 70db band, or any band for that matter. 
bigCount = 1;
sigCutoff = 0.05;
tarDB  = 3; %target db range, this I'm aiming for loudest sounds first.
toneTarget = 2; %period I wish to analyze.

slopeStore = [];
intStore = [];
slopeSpreadStore = [];
intSpreadStore = [];
valStore = [];
fullVals = [];
% counter = 1;
fullCount = 1;
fullType = [];
typeStore = [];

for i = 1:numFiles
    load(targetFiles{i})
    disp('Loading New File')
    disp(targetFiles{i})
    numUnits = length(s.DesignationName);
    for j = 1:numUnits
        disp(strcat('Examining Unit',s.DesignationName{j}))
        %now I need to go in and find significant responses that are also
        %positive
        fullType(fullCount) = masterData(j,7);
        
        sigVals = find(s.(s.DesignationName{j}).BinSigVals(2:end,tarDB,toneTarget) < sigCutoff);
        sigValsLaser = find(s.(s.DesignationName{j}).BinSigValsLaser(2:end,tarDB,toneTarget) < sigCutoff);
        posVals = find(s.(s.DesignationName{j}).BinDiff(2:end,tarDB,toneTarget) > 0);
        posValsLaser = find(s.(s.DesignationName{j}).BinDiffLaser(2:end,tarDB,toneTarget) > 0);
        %find intersects of both of these
        intersectNorm = intersect(sigVals,posVals);
        intersectLaser = intersect(sigValsLaser,posValsLaser);
        fullIntersect = intersect(intersectNorm,intersectLaser);
        disp(strcat('Found -',num2str(length(fullIntersect)),'Intersects'))
        normVals = s.(s.DesignationName{j}).BinDiff(2:end,tarDB,toneTarget);
        laserVals = s.(s.DesignationName{j}).BinDiffLaser(2:end,tarDB,toneTarget);
        fullVals{fullCount} = [normVals,laserVals];
        intersectStore(fullCount) = length(fullIntersect);
        fullCount = fullCount + 1;
        if length(fullIntersect) > 5
            disp('At least 5 significant responses in both laser and non laser')
            normVals = s.(s.DesignationName{j}).BinDiff(2:end,tarDB,toneTarget);
            laserVals = s.(s.DesignationName{j}).BinDiffLaser(2:end,tarDB,toneTarget);
            normVals = normVals(fullIntersect);
            laserVals = laserVals(fullIntersect);
            [b,bintr,bintjm] = gmregress(normVals,laserVals,sigCutoff);
            nameStore{bigCount} = strcat(targetFiles{i},s.DesignationName{j});
            typeStore(bigCount) = masterData(j,7);
            slopeStore(bigCount) = b(2);
            intStore(bigCount) = b(1);
            slopeSpreadStore(bigCount,:) = bintr(2,:);
            intSpreadStore(bigCount,:) = bintr(1,:);
            valStore{bigCount} = [normVals,laserVals];
            bigCount = bigCount + 1;
        end
    end
    
end

%we want to generate some plots of the overall dataset. Lets find an plot
%out MSNs and FSIs from the overall dataset for number of significant
%points.

msnSigNum = intersectStore(fullType==0);
pvSigNum = intersectStore(fullType==1);
msnSigNumPrune = msnSigNum(msnSigNum > 0);
pvSigNumPrune = pvSigNum(pvSigNum > 0);

hFig = figure;
set(hFig, 'Position', [80 80 500 800])
subplot(2,1,1)
hist(msnSigNumPrune,[1:1:16])
xlabel('Number of Significant Responses')
ylabel('Number of Units')
set(gca,'TickDir','out');
title(strcat('Sig Resp MSNs, Total:',num2str(length(msnSigNum)),'Non-Zero:',num2str(length(msnSigNumPrune))))

subplot(2,1,2)
hist(pvSigNumPrune,[1:1:16])
xlabel('Number of Significant Responses')
ylabel('Number of Units')
set(gca,'TickDir','out');
title(strcat('Sig Resp PVs, Total:',num2str(length(pvSigNum)),'Non-Zero:',num2str(length(pvSigNumPrune))))

spikeGraphName = '70DBMSNandPVNumberSigRespBothPos';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')


%find cell types
fsis = find(typeStore == 1);
msns = find(typeStore == 0);

%find significant slope changes. Trick here is to subtract by 1, which
%means that if the 95% confidence bounds include 1, there ends up being a
%negative on at least one side. Same applies for intercept, though dont
%need to subtract, since assumed intercept is zero. 

checkSigSlope = slopeSpreadStore - 1;
checkSigSlope = checkSigSlope(:,1) .*checkSigSlope(:,2);
checkSigSlope = find(checkSigSlope > 0);

checkSigInt = intSpreadStore(:,1) .* intSpreadStore(:,2);
checkSigInt = find(checkSigInt > 0);

%determine how many MSNs and fsis have significant changes
sigSlopeMSNs = intersect(msns,checkSigSlope);
sigSlopeFSIs = intersect(fsis,checkSigSlope);

sigIntMSNs = intersect(msns,checkSigInt);
sigIntFSIs = intersect(fsis,checkSigInt);

%now lets try plotting each line. The goal will be to make a big figure
%with multiple subplots. Each subplot will have black points, normal on the
%x, laser on the y, a unity line, and the regression line, in red. If
%regression slope is significant, or intercept is significant, this will be
%noted in the plot by the thickness of the line (both will produce a double
%width line). Units without sufficient responses will just be plotted as
%points. 


hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [80 80 1400 800])
for i = 1:length(msns)
    subplot(ceil(sqrt(length(msns))),ceil(sqrt(length(msns))),i)
    tarVals = valStore{msns(i)};
    minVal = min(min(tarVals));
    if minVal > 0
        minVal = 0;
    end
    maxVal = max(max(tarVals));
    plot(tarVals(:,1),tarVals(:,2),'k.')
    hold on
    plot([minVal maxVal],[minVal maxVal],'k')
    if ismember(msns(i),sigSlopeMSNs)
        plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',2)
    end
    if ismember(msns(i),sigIntMSNs)
        plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',1)
        plot(tarVals(:,1),tarVals(:,2),'r.')
    end
    plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',1)
%     title('')
    set(gca,'TickDir','out');
%     plot([minVal maxVal],[minVal maxVal],'k')
    axis equal
    xlim([minVal maxVal])
    ylim([minVal maxVal])
    
end
spikeGraphName = '70dbMSNIndivRegPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%and now for PVs

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [80 80 1400 800])
for i = 1:length(fsis)
    subplot(ceil(sqrt(length(fsis))),ceil(sqrt(length(fsis))),i)
    tarVals = valStore{fsis(i)};
    minVal = min(min(tarVals));
    if minVal > 0
        minVal = 0;
    end
    maxVal = max(max(tarVals));
    plot(tarVals(:,1),tarVals(:,2),'k.')
    hold on
    plot([minVal maxVal],[minVal maxVal],'k')
    if ismember(fsis(i),sigSlopeFSIs)
        plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',2)
    end
    if ismember(fsis(i),sigIntFSIs)
        plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',1)
        plot(tarVals(:,1),tarVals(:,2),'r.')
    end
    plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',1)
%     title('')
    set(gca,'TickDir','out');
%     plot([minVal maxVal],[minVal maxVal],'k')
    axis equal
    xlim([minVal maxVal])
    ylim([minVal maxVal])
    
end
spikeGraphName = '70dbPVIndivRegPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%% Now try for middle DB range

targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);


numFiles = length(targetFiles);

%okay, now that we have files, lets do something with them! 


%first thing i'm interested in checking is whether I can make significant
%regressions using the 70db band, or any band for that matter. 
bigCount = 1;
sigCutoff = 0.05;
tarDB  = 2; %target db range, this I'm aiming for loudest sounds first.
toneTarget = 2; %period I wish to analyze.

slopeStore = [];
intStore = [];
slopeSpreadStore = [];
intSpreadStore = [];
valStore = [];
fullVals = [];
% counter = 1;
fullCount = 1;
fullType = [];
typeStore = [];

for i = 1:numFiles
    load(targetFiles{i})
    disp('Loading New File')
    disp(targetFiles{i})
    numUnits = length(s.DesignationName);
    for j = 1:numUnits
        disp(strcat('Examining Unit',s.DesignationName{j}))
        %now I need to go in and find significant responses that are also
        %positive
        fullType(fullCount) = masterData(j,7);
        
        sigVals = find(s.(s.DesignationName{j}).BinSigVals(2:end,tarDB,toneTarget) < sigCutoff);
        sigValsLaser = find(s.(s.DesignationName{j}).BinSigValsLaser(2:end,tarDB,toneTarget) < sigCutoff);
        posVals = find(s.(s.DesignationName{j}).BinDiff(2:end,tarDB,toneTarget) > 0);
        posValsLaser = find(s.(s.DesignationName{j}).BinDiffLaser(2:end,tarDB,toneTarget) > 0);
        %find intersects of both of these
        intersectNorm = intersect(sigVals,posVals);
        intersectLaser = intersect(sigValsLaser,posValsLaser);
        fullIntersect = intersect(intersectNorm,intersectLaser);
        disp(strcat('Found -',num2str(length(fullIntersect)),'Intersects'))
        normVals = s.(s.DesignationName{j}).BinDiff(2:end,tarDB,toneTarget);
        laserVals = s.(s.DesignationName{j}).BinDiffLaser(2:end,tarDB,toneTarget);
        fullVals{fullCount} = [normVals,laserVals];
        intersectStore(fullCount) = length(fullIntersect);
        fullCount = fullCount + 1;
        if length(fullIntersect) > 5
            disp('At least 5 significant responses in both laser and non laser')
            normVals = s.(s.DesignationName{j}).BinDiff(2:end,tarDB,toneTarget);
            laserVals = s.(s.DesignationName{j}).BinDiffLaser(2:end,tarDB,toneTarget);
            normVals = normVals(fullIntersect);
            laserVals = laserVals(fullIntersect);
            [b,bintr,bintjm] = gmregress(normVals,laserVals,sigCutoff);
            nameStore{bigCount} = strcat(targetFiles{i},s.DesignationName{j});
            typeStore(bigCount) = masterData(j,7);
            slopeStore(bigCount) = b(2);
            intStore(bigCount) = b(1);
            slopeSpreadStore(bigCount,:) = bintr(2,:);
            intSpreadStore(bigCount,:) = bintr(1,:);
            valStore{bigCount} = [normVals,laserVals];
            bigCount = bigCount + 1;
        end
    end
    
end

%we want to generate some plots of the overall dataset. Lets find an plot
%out MSNs and FSIs from the overall dataset for number of significant
%points.

msnSigNum = intersectStore(fullType==0);
pvSigNum = intersectStore(fullType==1);
msnSigNumPrune = msnSigNum(msnSigNum > 0);
pvSigNumPrune = pvSigNum(pvSigNum > 0);

hFig = figure;
set(hFig, 'Position', [80 80 500 800])
subplot(2,1,1)
hist(msnSigNumPrune,[1:1:16])
xlabel('Number of Significant Responses')
ylabel('Number of Units')
set(gca,'TickDir','out');
title(strcat('Sig Resp MSNs, Total:',num2str(length(msnSigNum)),'Non-Zero:',num2str(length(msnSigNumPrune))))

subplot(2,1,2)
hist(pvSigNumPrune,[1:1:16])
xlabel('Number of Significant Responses')
ylabel('Number of Units')
set(gca,'TickDir','out');
title(strcat('Sig Resp PVs, Total:',num2str(length(pvSigNum)),'Non-Zero:',num2str(length(pvSigNumPrune))))

spikeGraphName = '60DBMSNandPVNumberSigRespBothPos';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')


%find cell types
fsis = find(typeStore == 1);
msns = find(typeStore == 0);

%find significant slope changes. Trick here is to subtract by 1, which
%means that if the 95% confidence bounds include 1, there ends up being a
%negative on at least one side. Same applies for intercept, though dont
%need to subtract, since assumed intercept is zero. 

checkSigSlope = slopeSpreadStore - 1;
checkSigSlope = checkSigSlope(:,1) .*checkSigSlope(:,2);
checkSigSlope = find(checkSigSlope > 0);

checkSigInt = intSpreadStore(:,1) .* intSpreadStore(:,2);
checkSigInt = find(checkSigInt > 0);

%determine how many MSNs and fsis have significant changes
sigSlopeMSNs = intersect(msns,checkSigSlope);
sigSlopeFSIs = intersect(fsis,checkSigSlope);

sigIntMSNs = intersect(msns,checkSigInt);
sigIntFSIs = intersect(fsis,checkSigInt);

%now lets try plotting each line. The goal will be to make a big figure
%with multiple subplots. Each subplot will have black points, normal on the
%x, laser on the y, a unity line, and the regression line, in red. If
%regression slope is significant, or intercept is significant, this will be
%noted in the plot by the thickness of the line (both will produce a double
%width line). Units without sufficient responses will just be plotted as
%points. 


hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [80 80 1400 800])
for i = 1:length(msns)
    subplot(ceil(sqrt(length(msns))),ceil(sqrt(length(msns))),i)
    tarVals = valStore{msns(i)};
    minVal = min(min(tarVals));
    if minVal > 0
        minVal = 0;
    end
    maxVal = max(max(tarVals));
    plot(tarVals(:,1),tarVals(:,2),'k.')
    hold on
    plot([minVal maxVal],[minVal maxVal],'k')
    if ismember(msns(i),sigSlopeMSNs)
        plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',2)
    end
    if ismember(msns(i),sigIntMSNs)
        plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',1)
        plot(tarVals(:,1),tarVals(:,2),'r.')
    end
    plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',1)
%     title('')
    set(gca,'TickDir','out');
%     plot([minVal maxVal],[minVal maxVal],'k')
    axis equal
    xlim([minVal maxVal])
    ylim([minVal maxVal])
    
end
spikeGraphName = '60dbMSNIndivRegPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%and now for PVs

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [80 80 1400 800])
for i = 1:length(fsis)
    subplot(ceil(sqrt(length(fsis))),ceil(sqrt(length(fsis))),i)
    tarVals = valStore{fsis(i)};
    minVal = min(min(tarVals));
    if minVal > 0
        minVal = 0;
    end
    maxVal = max(max(tarVals));
    plot(tarVals(:,1),tarVals(:,2),'k.')
    hold on
    plot([minVal maxVal],[minVal maxVal],'k')
    if ismember(fsis(i),sigSlopeFSIs)
        plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',2)
    end
    if ismember(fsis(i),sigIntFSIs)
        plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',1)
        plot(tarVals(:,1),tarVals(:,2),'r.')
    end
    plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',1)
%     title('')
    set(gca,'TickDir','out');
%     plot([minVal maxVal],[minVal maxVal],'k')
    axis equal
    xlim([minVal maxVal])
    ylim([minVal maxVal])
    
end
spikeGraphName = '60dbPVIndivRegPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%% Now try for low DB range

targets = what;
targetFiles = targets.mat;

masterIndex = strfind(targetFiles,'ML');
masterIndex = find(not(cellfun('isempty', masterIndex)));
targetFiles = targetFiles(masterIndex);


numFiles = length(targetFiles);

%okay, now that we have files, lets do something with them! 


%first thing i'm interested in checking is whether I can make significant
%regressions using the 70db band, or any band for that matter. 
bigCount = 1;
sigCutoff = 0.05;
tarDB  = 1; %target db range, this I'm aiming for loudest sounds first.
toneTarget = 2; %period I wish to analyze.

slopeStore = [];
intStore = [];
slopeSpreadStore = [];
intSpreadStore = [];
valStore = [];
fullVals = [];
% counter = 1;
fullCount = 1;
fullType = [];
typeStore = [];

for i = 1:numFiles
    load(targetFiles{i})
    disp('Loading New File')
    disp(targetFiles{i})
    numUnits = length(s.DesignationName);
    for j = 1:numUnits
        disp(strcat('Examining Unit',s.DesignationName{j}))
        %now I need to go in and find significant responses that are also
        %positive
        fullType(fullCount) = masterData(j,7);
        
        sigVals = find(s.(s.DesignationName{j}).BinSigVals(2:end,tarDB,toneTarget) < sigCutoff);
        sigValsLaser = find(s.(s.DesignationName{j}).BinSigValsLaser(2:end,tarDB,toneTarget) < sigCutoff);
        posVals = find(s.(s.DesignationName{j}).BinDiff(2:end,tarDB,toneTarget) > 0);
        posValsLaser = find(s.(s.DesignationName{j}).BinDiffLaser(2:end,tarDB,toneTarget) > 0);
        %find intersects of both of these
        intersectNorm = intersect(sigVals,posVals);
        intersectLaser = intersect(sigValsLaser,posValsLaser);
        fullIntersect = intersect(intersectNorm,intersectLaser);
        disp(strcat('Found -',num2str(length(fullIntersect)),'Intersects'))
        normVals = s.(s.DesignationName{j}).BinDiff(2:end,tarDB,toneTarget);
        laserVals = s.(s.DesignationName{j}).BinDiffLaser(2:end,tarDB,toneTarget);
        fullVals{fullCount} = [normVals,laserVals];
        intersectStore(fullCount) = length(fullIntersect);
        fullCount = fullCount + 1;
        if length(fullIntersect) > 5
            disp('At least 5 significant responses in both laser and non laser')
            normVals = s.(s.DesignationName{j}).BinDiff(2:end,tarDB,toneTarget);
            laserVals = s.(s.DesignationName{j}).BinDiffLaser(2:end,tarDB,toneTarget);
            normVals = normVals(fullIntersect);
            laserVals = laserVals(fullIntersect);
            [b,bintr,bintjm] = gmregress(normVals,laserVals,sigCutoff);
            nameStore{bigCount} = strcat(targetFiles{i},s.DesignationName{j});
            typeStore(bigCount) = masterData(j,7);
            slopeStore(bigCount) = b(2);
            intStore(bigCount) = b(1);
            slopeSpreadStore(bigCount,:) = bintr(2,:);
            intSpreadStore(bigCount,:) = bintr(1,:);
            valStore{bigCount} = [normVals,laserVals];
            bigCount = bigCount + 1;
        end
    end
    
end

%we want to generate some plots of the overall dataset. Lets find an plot
%out MSNs and FSIs from the overall dataset for number of significant
%points.

msnSigNum = intersectStore(fullType==0);
pvSigNum = intersectStore(fullType==1);
msnSigNumPrune = msnSigNum(msnSigNum > 0);
pvSigNumPrune = pvSigNum(pvSigNum > 0);

hFig = figure;
set(hFig, 'Position', [80 80 500 800])
subplot(2,1,1)
hist(msnSigNumPrune,[1:1:16])
xlabel('Number of Significant Responses')
ylabel('Number of Units')
set(gca,'TickDir','out');
title(strcat('Sig Resp MSNs, Total:',num2str(length(msnSigNum)),'Non-Zero:',num2str(length(msnSigNumPrune))))

subplot(2,1,2)
hist(pvSigNumPrune,[1:1:16])
xlabel('Number of Significant Responses')
ylabel('Number of Units')
set(gca,'TickDir','out');
title(strcat('Sig Resp PVs, Total:',num2str(length(pvSigNum)),'Non-Zero:',num2str(length(pvSigNumPrune))))

spikeGraphName = '50DBMSNandPVNumberSigRespBothPos';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')


%find cell types
fsis = find(typeStore == 1);
msns = find(typeStore == 0);

%find significant slope changes. Trick here is to subtract by 1, which
%means that if the 95% confidence bounds include 1, there ends up being a
%negative on at least one side. Same applies for intercept, though dont
%need to subtract, since assumed intercept is zero. 

checkSigSlope = slopeSpreadStore - 1;
checkSigSlope = checkSigSlope(:,1) .*checkSigSlope(:,2);
checkSigSlope = find(checkSigSlope > 0);

checkSigInt = intSpreadStore(:,1) .* intSpreadStore(:,2);
checkSigInt = find(checkSigInt > 0);

%determine how many MSNs and fsis have significant changes
sigSlopeMSNs = intersect(msns,checkSigSlope);
sigSlopeFSIs = intersect(fsis,checkSigSlope);

sigIntMSNs = intersect(msns,checkSigInt);
sigIntFSIs = intersect(fsis,checkSigInt);

%now lets try plotting each line. The goal will be to make a big figure
%with multiple subplots. Each subplot will have black points, normal on the
%x, laser on the y, a unity line, and the regression line, in red. If
%regression slope is significant, or intercept is significant, this will be
%noted in the plot by the thickness of the line (both will produce a double
%width line). Units without sufficient responses will just be plotted as
%points. 


hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [80 80 1400 800])
for i = 1:length(msns)
    subplot(ceil(sqrt(length(msns))),ceil(sqrt(length(msns))),i)
    tarVals = valStore{msns(i)};
    minVal = min(min(tarVals));
    if minVal > 0
        minVal = 0;
    end
    maxVal = max(max(tarVals));
    plot(tarVals(:,1),tarVals(:,2),'k.')
    hold on
    plot([minVal maxVal],[minVal maxVal],'k')
    if ismember(msns(i),sigSlopeMSNs)
        plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',2)
    end
    if ismember(msns(i),sigIntMSNs)
        plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',1)
        plot(tarVals(:,1),tarVals(:,2),'r.')
    end
    plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',1)
%     title('')
    set(gca,'TickDir','out');
%     plot([minVal maxVal],[minVal maxVal],'k')
    axis equal
    xlim([minVal maxVal])
    ylim([minVal maxVal])
    
end
spikeGraphName = '50dbMSNIndivRegPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%and now for PVs

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [80 80 1400 800])
for i = 1:length(fsis)
    subplot(ceil(sqrt(length(fsis))),ceil(sqrt(length(fsis))),i)
    tarVals = valStore{fsis(i)};
    minVal = min(min(tarVals));
    if minVal > 0
        minVal = 0;
    end
    maxVal = max(max(tarVals));
    plot(tarVals(:,1),tarVals(:,2),'k.')
    hold on
    plot([minVal maxVal],[minVal maxVal],'k')
    if ismember(fsis(i),sigSlopeFSIs)
        plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',2)
    end
    if ismember(fsis(i),sigIntFSIs)
        plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',1)
        plot(tarVals(:,1),tarVals(:,2),'r.')
    end
    plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',1)
%     title('')
    set(gca,'TickDir','out');
%     plot([minVal maxVal],[minVal maxVal],'k')
    axis equal
    xlim([minVal maxVal])
    ylim([minVal maxVal])
    
end
spikeGraphName = '50dbPVIndivRegPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%% now lets try this thing over again, with all data points (all freqs and dBs) but requiring matching still. 
bigCount = 1;
sigCutoff = 0.05;
tarDB  = 3; %target db range, this I'm aiming for loudest sounds first.
toneTarget = 2; %period I wish to analyze.
intersectStore = [];
slopeStore = [];
valStore = [];
intStore = [];
slopeSpreadStore = [];
intSpreadStore = [];
typeStore=[];

counter = 1;
for i = 1:numFiles
    load(targetFiles{i})
    disp('Loading New File')
    disp(targetFiles{i})
    numUnits = length(s.DesignationName);
    for j = 1:numUnits
        disp(strcat('Examining Unit',s.DesignationName{j}))
        %now I need to go in and find significant responses that are also
        %positive
        sigVals = find(s.(s.DesignationName{j}).BinSigVals(2:end,:,toneTarget) < sigCutoff);
        sigValsLaser = find(s.(s.DesignationName{j}).BinSigValsLaser(2:end,:,toneTarget) < sigCutoff);
        posVals = find(s.(s.DesignationName{j}).BinDiff(2:end,:,toneTarget) > 0);
        posValsLaser = find(s.(s.DesignationName{j}).BinDiffLaser(2:end,:,toneTarget) > 0);
        %find intersects of both of these
        intersectNorm = intersect(sigVals,posVals);
        intersectLaser = intersect(sigValsLaser,posValsLaser);
        fullIntersect = intersect(intersectNorm,intersectLaser);
        disp(strcat('Found -',num2str(length(fullIntersect)),'Intersects'))
        intersectStore(counter) = length(fullIntersect);
        counter = counter + 1;
        if length(fullIntersect) > 5
            disp('At least 5 significant responses in both laser and non laser')
            normVals = s.(s.DesignationName{j}).BinDiff(2:end,:,toneTarget);
            laserVals = s.(s.DesignationName{j}).BinDiffLaser(2:end,:,toneTarget);
            normVals = normVals(fullIntersect);
            laserVals = laserVals(fullIntersect);
            [b,bintr,bintjm] = gmregress(normVals,laserVals,sigCutoff);
            nameStore{bigCount} = strcat(targetFiles{i},s.DesignationName{j});
            typeStore(bigCount) = masterData(j,7);
            slopeStore(bigCount) = b(2);
            intStore(bigCount) = b(1);
            slopeSpreadStore(bigCount,:) = bintr(2,:);
            intSpreadStore(bigCount,:) = bintr(1,:);
            valStore{bigCount} = [reshape(normVals,[],1),reshape(laserVals,[],1)];
            bigCount = bigCount + 1;
        end
    end
    
end

%find cell types
fsis = find(typeStore == 1);
msns = find(typeStore == 0);

%find significant slope changes. Trick here is to subtract by 1, which
%means that if the 95% confidence bounds include 1, there ends up being a
%negative on at least one side. Same applies for intercept, though dont
%need to subtract, since assumed intercept is zero. 

checkSigSlope = slopeSpreadStore - 1;
checkSigSlope = checkSigSlope(:,1) .*checkSigSlope(:,2);
checkSigSlope = find(checkSigSlope > 0);

checkSigInt = intSpreadStore(:,1) .* intSpreadStore(:,2);
checkSigInt = find(checkSigInt > 0);

%determine how many MSNs and fsis have significant changes
sigSlopeMSNs = intersect(msns,checkSigSlope);
sigSlopeFSIs = intersect(fsis,checkSigSlope);

sigIntMSNs = intersect(msns,checkSigInt);
sigIntFSIs = intersect(fsis,checkSigInt);


hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [80 80 1400 800])
for i = 1:length(msns)
    subplot(ceil(sqrt(length(msns))),ceil(sqrt(length(msns))),i)
    tarVals = valStore{msns(i)};
    minVal = min(min(tarVals));
    if minVal > 0
        minVal = 0;
    end
    maxVal = max(max(tarVals));
    plot(tarVals(:,1),tarVals(:,2),'k.')
    hold on
    plot([minVal maxVal],[minVal maxVal],'k')
    if ismember(msns(i),sigSlopeMSNs)
        plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',2)
    end
    if ismember(msns(i),sigIntMSNs)
        plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',1)
        plot(tarVals(:,1),tarVals(:,2),'r.')
    end
    plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',1)
    
%     title('')
    set(gca,'TickDir','out');
%     plot([minVal maxVal],[minVal maxVal],'k')
    axis equal
    xlim([minVal maxVal])
    ylim([minVal maxVal])
    
end
spikeGraphName = 'alldbMSNIndivRegPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%and now for PVs

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [80 80 1400 800])
for i = 1:length(fsis)
    subplot(ceil(sqrt(length(fsis))),ceil(sqrt(length(fsis))),i)
    tarVals = valStore{fsis(i)};
    minVal = min(min(tarVals));
    if minVal > 0
        minVal = 0;
    end
    maxVal = max(max(tarVals));
    plot(tarVals(:,1),tarVals(:,2),'k.')
    hold on
    plot([minVal maxVal],[minVal maxVal],'k')
    if ismember(fsis(i),sigSlopeFSIs)
        plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',2)
    end
    if ismember(fsis(i),sigIntFSIs)
        plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',1)
        plot(tarVals(:,1),tarVals(:,2),'r.')
    end
    plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',1)
    
%     title('')
    set(gca,'TickDir','out');
%     plot([minVal maxVal],[minVal maxVal],'k')
    axis equal
    xlim([minVal maxVal])
    ylim([minVal maxVal])
    
end
spikeGraphName = 'alldbPVIndivRegPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')



%% now lets try this thing over again, with all data points (all freqs and
%dBs), without positive/intersect cutoff. 
bigCount = 1;
sigCutoff = 0.05;
tarDB  = 3; %target db range, this I'm aiming for loudest sounds first.
toneTarget = 2; %period I wish to analyze.
intersectStore = [];
slopeStore = [];
intStore = [];
slopeSpreadStore = [];
intSpreadStore = [];
bfStore = [];
valStore = [];
typeStore = [];

counter = 1;
for i = 1:numFiles
    load(targetFiles{i})
    disp('Loading New File')
    disp(targetFiles{i})
    numUnits = length(s.DesignationName);
    for j = 1:numUnits
        disp(strcat('Examining Unit',s.DesignationName{j}))
        %now I need to go in and find significant responses that are also
        %positive
        sigVals = find(s.(s.DesignationName{j}).BinSigVals(2:end,:,toneTarget) < sigCutoff);
        
        if length(sigVals) > 5
            disp('At least 5 significant responses overall')
            normVals = s.(s.DesignationName{j}).BinDiff(2:end,:,toneTarget);
            laserVals = s.(s.DesignationName{j}).BinDiffLaser(2:end,:,toneTarget);
            [b,bintr,bintjm] = gmregress(normVals,laserVals,sigCutoff);
            nameStore{bigCount} = strcat(targetFiles{i},s.DesignationName{j});
            typeStore(bigCount) = masterData(j,7);
            slopeStore(bigCount) = b(2);
            intStore(bigCount) = b(1);
            slopeSpreadStore(bigCount,:) = bintr(2,:);
            intSpreadStore(bigCount,:) = bintr(1,:);
            bfStore(bigCount,1) = masterData(j,12);
            bfStore(bigCount,2) = masterData(j,21);
            valStore{bigCount} = [reshape(normVals,[],1),reshape(laserVals,[],1)];
            bigCount = bigCount + 1;
        end
    end
    
end

%find cell types
fsis = find(typeStore == 1);
msns = find(typeStore == 0);

%find significant slope changes. Trick here is to subtract by 1, which
%means that if the 95% confidence bounds include 1, there ends up being a
%negative on at least one side. Same applies for intercept, though dont
%need to subtract, since assumed intercept is zero. 

checkSigSlope = slopeSpreadStore - 1;
checkSigSlope = checkSigSlope(:,1) .*checkSigSlope(:,2);
checkSigSlope = find(checkSigSlope > 0);

checkSigInt = intSpreadStore(:,1) .* intSpreadStore(:,2);
checkSigInt = find(checkSigInt > 0);

%determine how many MSNs and fsis have significant changes
sigSlopeMSNs = intersect(msns,checkSigSlope);
sigSlopeFSIs = intersect(fsis,checkSigSlope);

sigIntMSNs = intersect(msns,checkSigInt);
sigIntFSIs = intersect(fsis,checkSigInt);

%now lets find the ones that have differences in both?
doubleSigMSNs = intersect(sigSlopeMSNs,sigIntMSNs);
doubleSigFSIs = intersect(sigSlopeFSIs,sigIntFSIs);



hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [80 80 1400 800])
for i = 1:length(msns)
    subplot(ceil(sqrt(length(msns))),ceil(sqrt(length(msns))),i)
    tarVals = valStore{msns(i)};
    minVal = min(min(tarVals));
    if minVal > 0
        minVal = 0;
    end
    maxVal = max(max(tarVals));
    plot(tarVals(:,1),tarVals(:,2),'k.')
    hold on
    plot([minVal maxVal],[minVal maxVal],'k')
    if ismember(msns(i),sigSlopeMSNs) | ismember(msns(i),sigIntMSNs)
        plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',2)
    else
        plot([minVal maxVal],[minVal*slopeStore(msns(i))+intStore(msns(i)) maxVal*slopeStore(msns(i))+intStore(msns(i))],'r','LineWidth',1)
    end
%     title('')
    set(gca,'TickDir','out');
%     plot([minVal maxVal],[minVal maxVal],'k')
    axis equal
    xlim([minVal maxVal])
    ylim([minVal maxVal])
    
end
spikeGraphName = 'allMSNIndivRegPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')

%and now for PVs

hFig = figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.03 0.05], [0.03 0.01]);
set(hFig, 'Position', [80 80 1400 800])
for i = 1:length(fsis)
    subplot(ceil(sqrt(length(fsis))),ceil(sqrt(length(fsis))),i)
    tarVals = valStore{fsis(i)};
    minVal = min(min(tarVals));
    if minVal > 0
        minVal = 0;
    end
    maxVal = max(max(tarVals));
    plot(tarVals(:,1),tarVals(:,2),'k.')
    hold on
    plot([minVal maxVal],[minVal maxVal],'k')
    if ismember(fsis(i),sigSlopeFSIs) | ismember(fsis(i),sigIntFSIs)
        plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',2)
    else
        plot([minVal maxVal],[minVal*slopeStore(fsis(i))+intStore(fsis(i)) maxVal*slopeStore(fsis(i))+intStore(fsis(i))],'r','LineWidth',1)
    end
%     title('')
    set(gca,'TickDir','out');
%     plot([minVal maxVal],[minVal maxVal],'k')
    axis equal
    xlim([minVal maxVal])
    ylim([minVal maxVal])
    
end
spikeGraphName = 'allPVIndivRegPlot';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')





