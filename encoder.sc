%This is the sc file
%This is for reading the encoder data
%This has output A going to 7, B going to 8

int upA = 0
int vel = 0
int intWindow = 500

callback portin[7] up 
    upA = upA + 1
    if upA == 10 do
        vel = vel + 1
        portout[1] = 1
        do in 100
            portout[1] = 0
        end
        disp(vel)
        upA = 0
        do in intWindow
            portout[2] = 1
            do in 100
                portout[2] = 0
            end
            vel = vel - 1
            disp(vel)
        end
    end
end;
