clock(slave)

callback portin[7] up
    portout[5] = 1
    do in 500
        portout[5] = 0
    end
end;


