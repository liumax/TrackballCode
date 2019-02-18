function [tmf, xmf, sta_rtf, mid1_rtf, mid2_rtf] = calc_sta_mid_rtf(filter, filtopt)

if filtopt == 1
    gaussian = fspecial('gaussian', [2 2], 1);
end

% get parameters    
sprfile = filter.sprfile;
basestimfile = regexp(sprfile, '^\S+(?=(.spr))', 'match', 'once');
paramfile = [basestimfile '_param.mat'];
load(paramfile, 'MaxFM', 'MaxRD', 'faxis', 'taxis')

% get truncated faxis
faxis = filter.freq;

% get taxis
taxis = filter.time;

sta = filter.v_sta;
v1 = filter.v1;
v2 = filter.v2;

if filtopt == 1
    sta = imfilter(sta,gaussian);
    v1 = imfilter(v1,gaussian);
    v2 = imfilter(v2,gaussian);
end

[tmf, xmf, sta_rtf] = sta2rtf(sta, taxis, faxis, MaxFM, MaxRD);
[~, ~, mid1_rtf] = sta2rtf(v1, taxis, faxis, MaxFM, MaxRD);
[~, ~, mid2_rtf] = sta2rtf(v2, taxis, faxis, MaxFM, MaxRD);

end
    
    