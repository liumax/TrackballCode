tester = zeros(length(data),1);

sampRate = 192000;
targetTTL = 0.005; %target TTL duration in msec
whileCounter = 0;
lastPoint = 1;
while whileCounter == 0
    %generate random number to add to 300 ms
    randNum = round((rand(1)-0.5)*100);
    newInt = (500+randNum)/1000;
    %check to see if limit reached
    testEnd = lastPoint + newInt*sampRate + targetTTL*sampRate;
    if testEnd > length(tester)
        whileCounter = 1;
        break
    end
    %now add in ones to TTL
    tester(lastPoint+newInt*sampRate:lastPoint+newInt*sampRate+targetTTL*sampRate) = 1;
    lastPoint = lastPoint + newInt*sampRate;
    
end

data(:,2) = tester;