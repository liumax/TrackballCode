
tester = what;
tester = tester.mat;

fullCounter = 1;
bigStore = [];
for i = 1:length(tester)
    load(tester{i})
    numUnits = length(s.DesignationName);
    widthStore(fullCounter:fullCounter + numUnits - 1) = s.MasterSheet(:,2);
    modStore(fullCounter:fullCounter + numUnits - 1) = (s.MasterSheet(:,8)-s.MasterSheet(:,6))./(s.MasterSheet(:,6)+s.MasterSheet(:,8));
    distStore(fullCounter:fullCounter + numUnits - 1) = s.MasterSheet(:,4);
    fullCounter = fullCounter + numUnits;
end

%find MSNs and PVs
findMSNs = find(widthStore >= 0.0005);
findPVs = find(widthStore < 0.0004);

%adjust distance
newDist = distStore * -0.025 + 0.15;

dataset.findMSNs = findMSNs;
dataset.newDist = newDist;
dataset.modStore = modStore;

save('LaserDataset.mat','dataset')

figure
plot(dataset.newDist(dataset.findMSNs),dataset.modStore(dataset.findMSNs),'k.')
hold on
plot([0.15 1],[0 0],'b')
xlabel('Depth from Fiber Tip')
ylabel('Modulation Index')

%lets try and make averages across each site
xVal = dataset.newDist(dataset.findMSNs);
yVal = dataset.modStore(dataset.findMSNs);
numDepths = length(unique(xVal));
allDepths = unique(xVal);
meanStore = zeros(numDepths,1);
for i = 1:numDepths
    finder = find(xVal == allDepths(i));
    meanStore(i) = mean(yVal(finder));
end


%plot out distribution of units
figure
hist(xVal,[0.1:0.1:1])
xlim([0.1 1])
xlabel('Depth from Fiber Tip, mm')
ylabel('Number of Units')









