function [histog, bin_edges] = perievent_histog(spk_t, event_t, pre, pst, bin)
% Construct peri-event histograms
%
% Inputs
%	spk_t - vector of spike times 
%	event_t - vector of event times
%	pre - time before event at which to start histogram (negative value
%	assummed)
%	pst - time following event at which to stop histogram
%	bin - size of bins
%
% Output
%	histog - mean rate of spike events binned & expressed in spikes/sec
%	bin_edges - matching vector of edge time for histogram
%
% RST 05-08

bin_edges = pre:bin:pst;
hist_len = length(bin_edges);

n_events = length(event_t);
histog = zeros(1,hist_len);
bins = zeros(1,hist_len);
for i = 1:n_events
	histog = histog + histc( spk_t-event_t(i),bin_edges);
	% Correct n on bin-by-bin basis in case of early or late events (close
	% to edges
	bins = bins + ((bin_edges > -1*event_t(i))&(bin_edges < spk_t(end)-event_t(i)));
end
	
histog = histog ./ (bins .* bin);
histog(end) = [];
bin_edges(end) = [];
