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
plot([61+(i-1)*size(s.Processed.PhotoAverages,1) 61+(i-1)*size(s.Processed.PhotoAverages,1)],[minval maxval],'r')
end
