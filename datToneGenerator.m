clear all
close all
clc

[fname pname] = uiputfile('test1');

ao= analogoutput('nidaq','Dev2');
ch= addchannel(ao,0);
 

toneRange=[2000:2000:20000];
toneReps=20;
toneRecord = zeros(toneRange*toneReps,2);

%Should probably change to random presentation order!

t0 = clock;

disp(length(toneRange)*toneReps)

for i = 1:length(toneRange)
    disp(toneRange(i))
    for j = 1:toneReps
        sf = 100000;
        dur = 0.5;
        n = sf*dur;
        s=(1:n)/sf;
        s=sin(2*pi*toneRange(i)*s);
        sound(s,sf)
        toneRecord(j*i,1) = etime(clock,t0);
        toneRecord(j*1,2) = toneRange(i);
        putsample(ao,5)
        pause(0.1)
        putsample(ao,0)
        disp(i*j)
        pause(3)
    end
    pause(5)
end


save(fullfile(pname,fname),'toneRecord')

disp('ENDED')
