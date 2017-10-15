%this is code to run a three output pavlovian task. The goal here is to have catch trials, free trials, and normal cued trials for reward, nothing, and punishment. 

int itiDur
int soundPort %port for sound.
int outPort %port for reward/punishment
int outDur %duration for output
int toneOutDel %delay between tone and output

function 1 %This function serves to execute most of the code
    disp('TrialStart')
    do in itiDur
        %play sound
        disp('Tone Delivered')
        portout[soundPort] = 1
        portout[8] = 1
        do in 10
            portout[soundPort] = 0
            portout[8] = 1
        end
        do in toneOutDel
            portout[outPort] = 1
            do in outDur
                portout[outPort] = 0
            end
            do in 5000
                disp('PlotTime')
            end
        end
    end
end;

callback portin[2] up %lick has been detected
    disp('Lick Detected')
end;

