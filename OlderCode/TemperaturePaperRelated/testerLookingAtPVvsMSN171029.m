
%this is meant to look at some lumped data that has been combined on 10/29

msns = find(newSet(:,1) == 0);
pvs = find(newSet(:,1) == 1);

%can also define by reset boundaries.

msns = find(newSet(:,2) > 0.0005);
pvs = find(newSet(:,2) < 0.0004);


figure
hold on
plot([0:0.25:23],cumsum(pvRateHist)/sum(pvRateHist))
plot(mean(newSet(pvs,3)),0.5,'b.')
plot(median(newSet(pvs,3)),0.6,'b*')
plot([0:0.25:23],cumsum(msnRateHist)/sum(msnRateHist),'r')
plot(mean(newSet(msns,3)),0.5,'r.')
plot(median(newSet(msns,3)),0.6,'r*')
xlim([0 23])


