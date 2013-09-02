function events = GetNexStrobedEvents(nexname)
% function events = GetNexStrobedEvents(nexname)
%
% Get strobed events from NEX file
% (strobed events used for monkey data only)
% Parse events into "events" structure
%
%	Input- 
%		nexname - name of NEX file to read
%
%	Output-
%		events - structure containing times of parsed events
%
%	RST 2006-11-27



	[xa, xb, xc, ts, xd, char_vals] = nex_marker( nexname, 'Strobed', false);
	for s=1:length(char_vals)
		vals(s) = str2num(char_vals(s,:));
	end

	events = ParseStrobedEvents(ts,vals);
	
	