%this is code for reading out the TTLs from the Roland
%The goal is to detect all TTLs and respond to a subset of them.

clock(slave)

callback portin[2] up
	portout[1] = 1
	do in 10
		portout[1] = 0
	end
end;





