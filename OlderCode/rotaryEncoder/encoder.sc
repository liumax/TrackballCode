%This is the sc file
%This is for reading the encoder data
%This has output A going to 7, B going to 8

int counter = 0
int vel = 0
int intWindow = 500

callback portin[7] up 
    counter = counter + 1
    if counter == 3 do
        vel = vel + 1
        counter = 0
        do in intWindow
            vel = vel - 1
            if vel < 5 do %This is lower limit for running speed
                portout[1] = 0
            end
            if vel > 2 do
                portout[2] = 0
            end
        end
    end
    if vel > 6 do
        portout[1] = 1
    end
    if vel < 5 do
        portout[1] = 0
    end
    if vel < 1 do
        portout[2] = 1
    end
    if vel > 2 do
        portout[2] = 0
    end
end;
