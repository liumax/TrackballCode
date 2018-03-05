


%This uses the complete dataset that I generated!

%lets do some basic comparisons

%compare modulation of MSNs for significant differences

pvModARBR = (fullMasterARBR(fullMasterARBR(:,1)==1,8)-...
    fullMasterARBR(fullMasterARBR(:,1)==1,6))./...
    (fullMasterARBR(fullMasterARBR(:,1)==1,8)+...
    fullMasterARBR(fullMasterARBR(:,1)==1,6));

pvModIRES = (fullMasterIRES(fullMasterIRES(:,1)==1,8)-...
    fullMasterIRES(fullMasterIRES(:,1)==1,6))./...
    (fullMasterIRES(fullMasterIRES(:,1)==1,8)+...
    fullMasterIRES(fullMasterIRES(:,1)==1,6));

msnModARBR = (fullMasterARBR(fullMasterARBR(:,1)==0,8)-...
    fullMasterARBR(fullMasterARBR(:,1)==0,6))./...
    (fullMasterARBR(fullMasterARBR(:,1)==0,8)+...
    fullMasterARBR(fullMasterARBR(:,1)==0,6));

msnModIRES = (fullMasterIRES(fullMasterIRES(:,1)==0,8)-...
    fullMasterIRES(fullMasterIRES(:,1)==0,6))./...
    (fullMasterIRES(fullMasterIRES(:,1)==0,8)+...
    fullMasterIRES(fullMasterIRES(:,1)==0,6));

testPVMod = ranksum(pvModARBR,pvModIRES)
testMSNMod = ranksum(msnModARBR,msnModIRES)

%values are pvModTest 0.34406, msnModTest 8.36 e-6, so very significant

%now lets see if the distributions differ from noise

testPVNoiseARBR = signrank(pvModARBR)
testPVNoiseIRES = signrank(pvModIRES)
testMSNNoiseARBR = signrank(msnModARBR)
testMSNNoiseIRES = signrank(msnModIRES)

% testPVNoiseARBR =
% 
%     0.0345
% 
% 
% testPVNoiseIRES =
% 
%     0.6272
% 
% 
% testMSNNoiseARBR =
% 
%    3.2312e-05
% 
% 
% testMSNNoiseIRES =
% 
%    3.9806e-14

%so basically, the PV change in the ARBR is significantly different from a
%median of zero, not so for the IRES. Both MSN sets have medians
%significantly different from zero.


%look at restricted timing

pvResModARBR = (fullMasterARBR(fullMasterARBR(:,1)==1,11)-...
    fullMasterARBR(fullMasterARBR(:,1)==1,9))./...
    (fullMasterARBR(fullMasterARBR(:,1)==1,11)+...
    fullMasterARBR(fullMasterARBR(:,1)==1,9));

pvResModIRES = (fullMasterIRES(fullMasterIRES(:,1)==1,11)-...
    fullMasterIRES(fullMasterIRES(:,1)==1,9))./...
    (fullMasterIRES(fullMasterIRES(:,1)==1,11)+...
    fullMasterIRES(fullMasterIRES(:,1)==1,9));

msnResModARBR = (fullMasterARBR(fullMasterARBR(:,1)==0,11)-...
    fullMasterARBR(fullMasterARBR(:,1)==0,9))./...
    (fullMasterARBR(fullMasterARBR(:,1)==0,11)+...
    fullMasterARBR(fullMasterARBR(:,1)==0,9));

msnResModIRES = (fullMasterIRES(fullMasterIRES(:,1)==0,11)-...
    fullMasterIRES(fullMasterIRES(:,1)==0,9))./...
    (fullMasterIRES(fullMasterIRES(:,1)==0,11)+...
    fullMasterIRES(fullMasterIRES(:,1)==0,9));  


testPVResMod = ranksum(pvResModARBR,pvResModIRES)
testMSNResMod = ranksum(msnResModARBR,msnResModIRES)

% testPVResMod =
% 
%     0.4675
% 
% 
% testMSNResMod =
% 
%     0.0022
%so msn distributions are different by median, but PV not so.


%now lets see if the distributions differ from noise

testPVNoiseARBR = signrank(pvResModARBR)
testPVNoiseIRES = signrank(pvResModIRES)
testMSNNoiseARBR = signrank(msnResModARBR)
testMSNNoiseIRES = signrank(msnResModIRES)


% testPVNoiseARBR =
% 
%     0.0345
% 
% 
% testPVNoiseIRES =
% 
%     0.4576
% 
% 
% testMSNNoiseARBR =
% 
%     0.0162
% 
% 
% testMSNNoiseIRES =
% 
%    5.8057e-08

%again, PV Arbr has significantly different median, not so for IRES.
%for MSNs, looks like both are significantly different. 

%figure out percent suppression of PV cells in both cases. Use modulation
%of -0.3 as cutoff

perModPVARBR = length(find(pvModARBR<-0.3))/length(pvModARBR)
perModPVIRES = length(find(pvModIRES<-0.3))/length(pvModARBR)

% perModPVARBR =
% 
%     0.2963
% 
% 
% perModPVIRES =
% 
%     0.2222
%of course, this difference appears not to be significantly different by
%median...

%what if we select a harsher line? 
perModPVARBR = length(find(pvModARBR<-0.5))/length(pvModARBR)
perModPVIRES = length(find(pvModIRES<-0.5))/length(pvModARBR)

% perModPVARBR =
% 
%     0.1852
% 
% 
% perModPVIRES =
% 
%     0.2222
%equalizes a bit, brings down proportion of ARBR cells modulated.

%okay what about firing rate distributions?


meanrateMSNARBR = fullMasterARBR(fullMasterARBR(:,1)==0,6);
meanrateMSNIRES = fullMasterIRES(fullMasterIRES(:,1)==0,6);
meanratePVARBR = fullMasterARBR(fullMasterARBR(:,1)==1,6);
meanratePVIRES = fullMasterIRES(fullMasterIRES(:,1)==1,6);

%stats tests!

testMeanRateMSN = ranksum(meanrateMSNARBR,meanrateMSNIRES)
testMeanRatePV = ranksum(meanratePVARBR,meanratePVIRES)

% testMeanRateMSN =
% 
%     0.1611
% 
% 
% testMeanRatePV =
% 
%     0.0686
%no significant differences...

%try ks stats test instead
% 
% testKSMeanRateMSN = kstest(meanrateMSNARBR,meanrateMSNIRES)
% testKSMeanRatePV = kstest(meanratePVARBR,meanratePVIRES)

%derp derp derp, doesnt work! XD moving on now.....


%how about we plot vs distance?
figure
plot(fullMasterARBR(fullMasterARBR(:,1)==0,4),msnModARBR,'k.')
hold on
plot(fullMasterIRES(fullMasterIRES(:,1)==0,4),msnModIRES,'r.')

%linear regression
tester = [ones(length(msnModARBR),1),(fullMasterARBR(fullMasterARBR(:,1)==0,4)/-1000)]\msnModARBR;
testCalc = [ones(length(msnModARBR),1),(fullMasterARBR(fullMasterARBR(:,1)==0,4)/-1000)]*tester;
figure
plot(fullMasterARBR(fullMasterARBR(:,1)==0,4)/-1000,msnModARBR,'k.')
hold on
plot(fullMasterARBR(fullMasterARBR(:,1)==0,4)/-1000,testCalc,'r')
set (gca,'Xdir','reverse')

Rsq1 = 1 - sum((msnModARBR - testCalc).^2)/sum((msnModARBR - mean(msnModARBR)).^2);

