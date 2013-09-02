function [sdf, bin_edges] = perievent_sdf(spk_t, event_t, pre, pst, smooth)
% Construct peri-event SDF from vector of spike times
%
% Inputs
%	spk_t - vector of spike times
%	event_t - vector of event times
%	pre - time before event at which to start histogram (negative value
%	assummed)
%	pst - time following event at which to stop histogram
%	smooth - sigma of gaussian (in msec)
%
%
% Output
%	sdf - mean rate of spike events convolved w/ gaussian kernel & expressed in spikes/sec
%	bin_edges - matching vector of independent sample intervals
%
% RST 2005-08

bin_edges = pre:0.001:pst;
data_len = length(bin_edges);

raster = perievent_raster(spk_t, event_t, pre, pst);
delta = spk_t2delta( raster, data_len);
sdf = delta2sdf(delta,smooth);

