%port 1 is initial trigger from noldus
%port 2 is noldus syncing
%port 3 is output to optoscope control (trigger on) AND input for sync
%port 4 is output to amber LED

int itiTime = 60000




function 1
    do in itiTime
        portout[4] = 1
        do in itiTime
            portout[4] = 0
            trigger(1)
        end
    end
end;

callback portin[1] up
    portout[3] = 1
    trigger(1)
end;

