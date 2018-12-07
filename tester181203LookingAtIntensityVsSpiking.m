
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
            [Y,I] = max(s.(s.DesignationName{i}).BinTone(:,9)); %I is best frequency
            if I
                intensityResponse(:,i) = s.(s.DesignationName{i}).BinTone(I,:);
                intensityResponseLaser(:,i) = s.(s.DesignationName{i}).BinToneLaser(I,:);
                bigInt(:,counter) = s.(s.DesignationName{i}).BinTone(I,:);
                bigIntLaser(:,counter) = s.(s.DesignationName{i}).BinToneLaser(I,:);
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
