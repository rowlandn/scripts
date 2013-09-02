function [rstx, rsty] = spk_t2rastervecs(spk_t, y, tickheight)
% Make x & y vectors for raster plotting from spike times
%
%function [rstx, rsty] = spk_t2rastervecs(spk_t,tickheight)
%
% Inputs:
%		spk_t - times of spikes
%		y - Desired y-position of raster
%		tickheight - y-vector is made this height
len = length(spk_t);

rstx = zeros(2,len);
rstx(1,:) = spk_t;
rstx(2,:) = spk_t;
rsty = zeros(2,len);
rsty(1,:) = y;
rsty(2,:) = y+tickheight;


