tester = s.Processed.PhotoAverages(:,:,end);
testLin = reshape(tester,1,[]);
newLin = zeros(length(testLin),1);
for i = 1:size(s.Processed.PhotoAverages,2)
    newLin(size(s.Processed.PhotoAverages,1)*(i-1)+1:size(s.Processed.PhotoAverages,1)*(i)) = tester(:,i) - tester(1,i);
end
figure
plot(testLin)

figure
hold on
plot(newLin)
maxval = max(max(testLin));
minval = min(min(testLin));
for i = 1:size(s.Processed.PhotoAverages,2)
plot([i*size(s.Processed.PhotoAverages,1) i*size(s.Processed.PhotoAverages,1)],[minval maxval],'k')
plot([143+(i-1)*size(s.Processed.PhotoAverages,1) 143+(i-1)*size(s.Processed.PhotoAverages,1)],[minval maxval],'r')
end

%for inscopix
tester = s.Processed.PhotoAverages;
testLin = reshape(tester,1,[]);
newLin = zeros(length(testLin),1);
for i = 1:size(s.Processed.PhotoAverages,2)
    newLin(size(s.Processed.PhotoAverages,1)*(i-1)+1:size(s.Processed.PhotoAverages,1)*(i)) = tester(:,i) - mean(tester(50:59,i));
end
% figure
% plot(testLin)

windowSize = -s.Parameters.RasterWindow(1)/(s.Parameters.RasterWindow(2) - s.Parameters.RasterWindow(1));
windowSize = windowSize * size(s.Processed.PhotoAverages,1);
hFig = figure
set(hFig, 'Position', [10 80 1240 850])
hold on
plot(newLin)
maxval = max(max(testLin));
minval = min(min(testLin));
for i = 1:size(s.Processed.PhotoAverages,2)
plot([i*size(s.Processed.PhotoAverages,1) i*size(s.Processed.PhotoAverages,1)],[minval maxval],'k','LineWidth',2)
plot([windowSize+(i-1)*size(s.Processed.PhotoAverages,1) windowSize+(i-1)*size(s.Processed.PhotoAverages,1)],[minval maxval],'r')
end
xlim([1 length(newLin)])
ylim([minval maxval])
ylabel('Frequencies')
xlabel('Fluorescence')
title('Inscopix Photometry Tuning Curve')
savefig(hFig,'inscopixPhotometryTuningLinear');

%save as PDF with correct name
set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hFig,'inscopixPhotometryTuningLinear','-dpdf','-r0')


