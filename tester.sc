%These are variables for tracking running disk.
int rando = 0
int intWindow = 500
int velo = 0


callback portin[7] up 
    rando = rando + 1
    if rando == 3 do
        velo = velo + 1
        disp(velo)
        rando = 0
        do in intWindow
            velo = velo - 1 %If I eliminate this, the code works fine.
            disp(velo)
        end
    end
end;