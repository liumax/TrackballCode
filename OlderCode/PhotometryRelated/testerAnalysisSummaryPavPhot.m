%This is code meant to analyze the outputs of analysisWrapperPavPhot


tester = what;
fileNames = tester.mat;


for i = 1:length(fileNames)
    %load file
    load(fileNames{i})
    bigResp(i,:) = keyVals(:,2);
    smallResp(i,:) = keyVals(:,3);
    
end


subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
hFig = figure
set(hFig, 'Position', [10 80 1240 850])
subplot(2,3,1)
plot(bigResp')
title('Big Response Changes')
subplot(2,3,4)
plot(bigResp(:,[1,3])')
title('Big Response Changes Beginning and End')
subplot(2,3,2)
plot(smallResp')
title('Small Response Changes')

subplot(2,3,5)
plot(smallResp(:,[1,3])')
title('Small Response Changes Beginning and End')
subplot(2,3,3)
plot((bigResp./smallResp)')
title('Ratio, Big/Small, all sessions')
subplot(2,3,6)
plot((bigResp(:,[1,3])./smallResp(:,[1,3]))')
title('Ratio, Big/Small, all sessions')







