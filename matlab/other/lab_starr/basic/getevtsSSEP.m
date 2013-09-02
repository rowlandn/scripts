function [events,values] = getevtsSSEP(time, chan, rng1, rng2)
% function [events,values] = getevts(time, chan, inds1, inds2)
%
%	Function to find times at which chan steps from
%		rng1 to rng2
%
% RST 2006-04-05
% SS 10-16-2008: edited to eliminate TIME_RNG and only find inds where
% voltage is within rng2.  Used to perform SSEP analysis

TIME_RNG = 0;		% msec allowed for transition
initial_event = false;


if ~isnan(rng1)
	in1 = inrange(chan(1:end-TIME_RNG),rng1);
	in2 = inrange(chan(1+TIME_RNG:end),rng2);
else
	in1 = outrange(chan(1:end-TIME_RNG),rng2);
	in2 = inrange(chan(1+TIME_RNG:end),rng2);
	
	if inrange(chan(1),rng2)
		initial_event = true;
	end
end

% inds = find( in1 & in2 );
inds = find( in2 );
last_inds = find( diff(inds)>1 );

if ~isempty(inds)
	if initial_event
		events = time( [1 inds(last_inds) inds(end)] );
		values = chan( [1 inds(last_inds) inds(end)] );
	else
		events = time( [inds(last_inds)' inds(end)'] );
		values = chan( [inds(last_inds)' inds(end)'] );
	end
else
		events = [];
		values = [];
end

