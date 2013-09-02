function raster = perievent_raster(spk_t, event_t, pre, pst)
% PERI-EVENT RASTER
% function raster = perievent_raster(spk_t, event_t, pre, pst)
% 
%	Make peri-event raster array of spike times - 1 row per trial
%	NaN pads variable length rows

n_events = length(event_t);
raster = NaN*zeros(1,1);
rst_len = size(raster,2);
for i = 1:n_events
	pret = event_t(i)+pre;
	pstt = event_t(i)+pst;
	rstt = spk_t( find( spk_t>=pret & spk_t<pstt ) ) - event_t(i);
	rsttlen = size(rstt,2);
	% Resize raster array if necessary
	if rsttlen > rst_len
		raster(:,(rst_len+1:rsttlen)) = NaN;
		rst_len = size(raster,2);
	elseif rsttlen < rst_len
		rstt(rsttlen+1:rst_len) = NaN;
	end
	raster = [ raster; rstt ];
end
% Delete 1st dummy row when finished
raster(1,:) = [];
