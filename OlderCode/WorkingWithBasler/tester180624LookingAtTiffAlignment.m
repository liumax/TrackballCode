fileName = '180623_ML180515F_IRsyncTest';
extractTimeBinaryFile(fileName);
extractDioBinaryFiles(fileName);


[DIOData] = readTrodesExtractedDataFile('180623_ML180515F_IRsyncTest.dio_MCU_Aux1.dat');
sampleRate = 30000;
dioTimes = double(DIOData.fields(1).data)/sampleRate;
dioState = double(DIOData.fields(2).data);
frameTimes = [dioTimes,dioState];

[DIOData] = readTrodesExtractedDataFile('180623_ML180515F_IRsyncTest.dio_MCU_Din5.dat');
dioTimes = double(DIOData.fields(1).data)/sampleRate;
%extracts states. States should be 0 or 1.
dioState = double(DIOData.fields(2).data);
lightsOn = find(dioState == 1);

lightsOnTime = dioTimes(lightsOn);

matchTimes = zeros(200,1);
for i = 1:200
    matchTimes(i) = find(trueFrames(:,1) < lightsOnTime(i),1,'last');
end
%curate. If frame time is more than 2 msec before light onset, default to
%next frame
findBig = find(trueFrames(matchTimes)>0.002);
matchTimes(findBig) = matchTimes(findBig) + 1;

tester = dir('*.tiff');
tiffNames = {tester.name};
[cs,index] = sort_nat(tiffNames);

for i = 1:length(matchTimes)
    
    BIGctrlIm = imread(cs{matchTimes(i)-2})-imread(cs{matchTimes(i)-3});
    BIGctrlImStore{i} = BIGctrlIm(4:20,110:130,:);
    BIGctrlMean(i) = mean(mean(mean(BIGctrlIm(4:20,110:130,:))));
    
    ctrlIm = imread(cs{matchTimes(i)-1})-imread(cs{matchTimes(i)-2});
    ctrlImStore{i} = ctrlIm(4:20,110:130,:);
    ctrlMean(i) = mean(mean(mean(ctrlIm(4:20,110:130,:))));
    testIm = imread(cs{matchTimes(i)})-imread(cs{matchTimes(i)-1});
    testImStore{i} = testIm(4:20,110:130,:);
    testMean(i) = mean(mean(mean(testIm(4:20,110:130,:))));
    
    futestIm = imread(cs{matchTimes(i)+1})-imread(cs{matchTimes(i)});
    futestImStore{i} = futestIm(4:20,110:130,:);
    futestMean(i) = mean(mean(mean(futestIm(4:20,110:130,:))));
end


subplot = @(m,n,p) subtightplot (m, n, p, [0.01], [0.01 0.01], [0.01 0.01]);


figure
for i = 1:length(matchTimes)
    subplot(10,20,i)
    imagesc(BIGctrlImStore{i})
%     title(num2str(i))
end


% subplot = @(m,n,p) subtightplot (m, n, p, [0.01], [0.01 0.01], [0.01 0.01]);
figure
for i = 1:length(matchTimes)
    subplot(10,20,i)
    imagesc(ctrlImStore{i})
%     title(num2str(i))
end



figure
for i = 1:length(matchTimes)
    subplot(10,20,i)
    imagesc(testImStore{i})
%     title(num2str(i))
end



figure
for i = 1:length(matchTimes)
    subplot(10,20,i)
    imagesc(futestImStore{i})
%     title(num2str(i))
end



figure
plot(testMean)
hold on
plot(futestMean,'r')
plot(ctrlMean,'g')


timeDiffs = lightsOnTime - trueFrames(matchTimes,1);

findSmall = find(timeDiffs < 0.003);


figure
plot(testMean)
hold on
plot(findSmall,testMean(findSmall),'c*')
plot(futestMean,'r')
plot(ctrlMean,'g')
plot(testMean + futestMean,'m')

