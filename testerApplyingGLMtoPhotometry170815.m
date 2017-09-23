%%This code is meant to look at an analysis file and perform glmfit to the
%%data. 


rasterWindow = s.Parameters.RasterWindow;

newWindow = [-1,3];
timeBin = 0.1;

%pull photometry traces
toneRaster = s.PhotoRaster.ToneRaster;
rasterTimes = zeros(length(toneRaster),1);
rasterTimes(1) = rasterWindow(1);
for i=2:length(toneRaster)
    rasterTimes(i) = rasterTimes(i-1) + s.PhotoRaster.TimeStep;
end

%convert to new timescale
newVector = [newWindow(1):timeBin:newWindow(2)];

for i = 1:size(toneRaster,2)
    newPhot(:,i) = interp1(rasterTimes,toneRaster(:,i),newVector);
end

%pull velocity traces
velRaster = s.VelRaster.ToneRaster;
velTimes = s.VelRaster.Axis;


for i = 1:size(velRaster,2)
    newVel(:,i) = interp1(velTimes,velRaster(:,i),newVector);
end

%now I need to pull tone order
toneOrder = zeros(length(newVector),length(s.MBED.ToneDelivery));

toneOrder(:,s.MBED.HiTrials) = 1;

%now pull licking data

lickRasters = s.Licking.ToneRaster;
%preprocess by removing excess values
lickRasters(lickRasters(:,1) < newWindow(1),:) = [];
lickRasters(lickRasters(:,1) > newWindow(2),:) = [];

for i = 1:length(s.MBED.ToneDelivery)
    findLicks = find(lickRasters(:,2) == i);
    if length(findLicks) >0
        newLicks(:,i) = hist(lickRasters(findLicks,1),newVector);
    else
        newLicks(:,i) = zeros(length(newVector),1);
    end
end


%now perform GLM fits

for i = 1:length(newVector)
    X(:,1) = toneOrder(i,:);
    X(:,2) = newVel(i,:);
    X(:,3) = newLicks(i,:);
    X = zscore(X);
    %calculate variance inflation factor
    R0 = corrcoef(X); % correlation matrix
    V=diag(inv(R0))';
%     Variance inflation factor (VIF) quantifies how much the variance is 
% inflated due to collinearity of regressor matrix columns. 
% i_th entry in the output vector is the variance inflation factor 
% for the i_th predictor, which indicates how much the variance of 
% the i_th predictor is inflated due to collinearity.
    %store
    varStore(:,i) = V;
    y = newPhot(i,:);
    y = zscore(y);
    [beta non stats] = glmfit(X,y);
    betaStore(:,i) = beta;
    statStore(:,i) = stats.p;
end

hFig= figure
subplot(2,1,1)
hold on
plot(betaStore(2:end,:)')
legend('ToneOrder','Velocity','Licks')
for i = 2:size(betaStore,1)
    findSig = find(statStore(i,:) < 0.05);
    plot(findSig,betaStore(i,findSig),'r*')
end
set(gca,'XTick',[1:10:length(newVector)]);
set(gca,'XTickLabel',[newWindow(1):newWindow(2)]);
title(strcat('GLM Fit for',testNames{bigInd}))
xlabel('Time(s)')
ylabel('Beta Coefficient')



subplot(2,1,2)
plot(varStore')
set(gca,'XTick',[1:10:length(newVector)]);
set(gca,'XTickLabel',[newWindow(1):newWindow(2)]);
title('Variance Inflation Factor')
xlabel('Time(s)')
ylabel('VIF')

spikeGraphName = strcat(testNames{bigInd}(1:end-4),'GLMFigure');
savefig(hFig,spikeGraphName);

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,spikeGraphName,'-dpdf','-r0')


X = [];
varStore = [];
betaStore = [];
newVel = [];
newPhot = [];
newLicks = [];