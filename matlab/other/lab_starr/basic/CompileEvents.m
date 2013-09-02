function trials = CompileEvents(home_ev,light_ev,move_ev,...
	other_events,target)
%
%
%

if isempty(light_ev)
	trials.home = [];
	trials.light = [];
	trials.move = [];
	trials.target = [];
end

for i=1:length(light_ev)
	
	other_evs = [other_events, move_ev];
	% Find home event immediately preceeding right light
	hm_ev = home_ev( max( find( home_ev<light_ev(i) ) ) );
	
	% Make sure no other events intervene between home & right event
	er_ind = find( other_evs>hm_ev & other_evs<light_ev(i) );
	if ~isempty(er_ind)
		error('cannot parse trials correctly');
	end
	
	other_evs = [other_events, home_ev];
	% Find move event immediately following light
	mv_ev = move_ev( min( find( move_ev>light_ev(i) ) ) );
	
	if ~isempty(mv_ev)
		% Make sure no other events intervene between home & right event
		er_ind = find( other_evs>light_ev(i) & other_evs<mv_ev );
		if ~isempty(er_ind)
			error('cannot parse trials correctly');
		end
	else
		mv_ev = NaN;
	end
	
	trials.home(i) = hm_ev;
	trials.light(i) = light_ev(i);
	trials.move(i) = mv_ev;
	trials.target(i) = target;
end

	
	