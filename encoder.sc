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
        disp(vel)
        upA = 0
        do in intWindow
            vel = vel - 1
            disp(vel)
        end
    end
end;
