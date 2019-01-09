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
counter = 1;
fullCount = 1;
fullType = [];

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
        fullCount = fullCount + 1;
        sigVals = find(s.(s.DesignationName{j}).BinSigVals(2:end,tarDB,toneTarget) < sigCutoff);
        sigValsLaser = find(s.(s.DesignationName{j}).BinSigValsLaser(2:end,tarDB,toneTarget) < sigCutoff);
        posVals = find(s.(s.DesignationName{j}).BinDiff(2:end,tarDB,toneTarget) > 0);
        posValsLaser = find(s.(s.DesignationName{j}).BinDiffLaser(2:end,tarDB,toneTarget) > 0);
        %find intersects of both of these
        intersectNorm = intersect(sigVals,posVals);
        intersectLaser = intersect(sigValsLaser,posValsLaser);
        fullIntersect = intersect(intersectNorm,intersectLaser);
        disp(strcat('Found -',num2str(length(fullIntersect)),'Intersects'))
        intersectStore(counter) = length(fullIntersect);
        counter = counter + 1;
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



%now lets try this thing over again, with all data points (all freqs and dBs) 
bigCount = 1;
sigCutoff = 0.05;
tarDB  = 3; %target db range, this I'm aiming for loudest sounds first.
toneTarget = 2; %period I wish to analyze.
intersectStore = [];
slopeStore = [];
intStore = [];
slopeSpreadStore = [];
intSpreadStore = [];

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




%now lets try this thing over again, with all data points (all freqs and
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
            valStore(bigCount,1) = mean(mean(normVals));
            valStore(bigCount,2) = mean(mean(laserVals));
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

%lets try and categorize these units?

hFig = figure;
plot(valStore(msns,1),valStore(msns,2),'r.')
hold on
plot([-0.5 max(max(valStore(msns)))],[-0.5 max(max(valStore(msns)))],'b')
xlim([-0.5 2])
ylim([-0.5 2])
xlabel('No Laser Average')
ylabel('Laser Average')

spikeGraphName = 'msnSpikeScatterAverageAllResponses';
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-djpeg','-r0')









