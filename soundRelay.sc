
clock(slave)

callback portin[2] up
    portout[1] = 1
    do in 100
        portout[1] = 0
    end
end;



