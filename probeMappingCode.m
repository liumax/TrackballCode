channelNum = 32;

%generates empty array
orderVector = zeros(channelNum,1);

%places neuroTech designations in the correct order. Here, you want shank
%designation to be organized into the actual physical orientation of the
%shanks, numbered 1-32
for i = 1:channelNum
    x = find(shankdesigByPosition == i);
    orderVector(i) = shankdesigByNeuroTech(x);
end
%generates empty array for adaptor map. Here, you want adaptor mapped as
%though you are looking directly at the adaptor. R, G, and blanks, are
%labeled as zeros.
adaptorVector = zeros(channelNum,1);
%finds correct values in neurotech array, maps onto neuronexus array
for i = 1:channelNum
    x = find(adaptorByNeurotech == orderVector(i));
    adaptorVector(i) = adaptorByNexus(x);
end

%generates empty array for omnetics map. Have this be an 2x16 array without
%the R or G, since those are unnecessary. This should output values between
%0 and 31.
omneticsVector = zeros(channelNum,1);
%maps these values into trodes numbers
for i = 1:channelNum
    x = find(omneticsByNexus == adaptorVector(i))
    omneticsVector(i) = omneticsByTrodes(i);
end
