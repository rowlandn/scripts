function events = ParseStrobedEvents(ts,val)
% function events = ParseStrobedEvents(ts,val)
%
%	Written specifically for dystonia task
%	receives timestampes (ts) and strobe values 
%	returns structure of task events broken into trials
%
%	RST 2006-01-09

VERBOSE=false;

% Definitions of dystonia task codes
FIX	= 1;
CENTER_TOUCH = 2;	%2
RIGHT_TARGET = 4;	%3;
LEFT_TARGET = 8;	%4;
CENTER_LEAVE = 16;	%5;
RIGHT_TOUCH = 32;	%6;
LEFT_TOUCH = 64;	%7;
% Skip bits 8 & 9 'case they're used for spks in TAP files
REWARD = 512;		%10;
BEGIN = 1024;			% 11;
END = 2048;			% 12;

% Bit numbers
EV_B = 13;		% 4096 Bit denoting event codes (versus codes for trial # & start/stop)
TR_END = 14;	% 20000 - End of trial
TR_SAV = 15;	% 40000 - trial saved by HEBB

LEFT = 0;
RIGHT = 1;

evt_inds = find( bitget(val,EV_B) );
start_inds = find( ~bitget(val,EV_B) & ~bitget(val,TR_END));
stop_inds = find( ~bitget(val,EV_B) & bitget(val,TR_END));
nval = length(val);

% Find matching trial starts and stops
for i=1:length(start_inds)
	si = start_inds(i);
	start_i(i) = si;
	start_t(i) = ts(si);
	stops = find( ts(stop_inds)>start_t(i) );	%stops always follow starts
	if ~isempty(stops)
		stop_i(i) = min( stop_inds(stops));
		stop_t(i) = ts( stop_i(i) );	% assign stop time		
	else
		start_i(i) = [];	% else delete last start if no matching stop found
		start_t(i) = [];
	end
end
% check to make sure starts & stops make sense
err = ( length(start_i) ~= length(start_i) );
err = err | ( length(start_t) ~= length(start_t) );
if ~err
	for i=1:length(start_i)
		err = err | ( stop_t(i) < start_t(i) );
	end
end
if err
	error(['Start/Stop timing mismatch discovered in dysEvent2Trial']);
end

events.ntrials = length(start_i);
nulra = NaN * zeros(1,events.ntrials);
events.fix=nulra;
events.home=nulra;
events.light=nulra;
events.target= 46 + zeros(1,events.ntrials);	% '.' = ascii 46 - used to mark missing val
events.move=nulra;
events.touch=nulra;
events.touchDir= 46 + zeros(1,events.ntrials);	% '.' = ascii 46 - used to mark missing val
events.reward=nulra;
events.success= false(1,events.ntrials);

% strip off event bits
evts = bitset( val, EV_B, 0);

% Now populate trials with task event information
for t=1:events.ntrials
	% get event indices for current trial
	tr_i = evt_inds( find(evt_inds>start_i(t) & evt_inds<stop_i(t)) );	
	
	events.start(t) = start_t(t);
	events.stop(t) = stop_t(t);
	events.taptrial(t) = val( start_i(t) );
	for i=1:length(tr_i)
		ind = tr_i(i);
		switch evts( ind )
			case FIX
				events.fix(t) = ts(ind);
			case CENTER_TOUCH
				events.home(t) = ts(ind);
			case LEFT_TARGET
				events.light(t) = ts(ind);
				events.target(t) = 'L';
			case RIGHT_TARGET
				events.light(t) = ts(ind);
				events.target(t) = 'R';
			case CENTER_LEAVE
				events.move(t) = ts(ind);
			case LEFT_TOUCH
				events.touch(t) = ts(ind);
				events.touchDir(t) = 'L';
			case RIGHT_TOUCH
				events.touch(t) = ts(ind);
				events.touchDir(t) = 'R';
			case REWARD
				events.reward(t) = ts(ind);
				if events.touch(t) & ( events.touchDir(t) == events.target(t) )
					events.success(t) = true;
				end
			case BEGIN
				events.beg(t) = ts(ind);
			case END
				events.end(t) = ts(ind);
			otherwise
				warning(['Detected an unknown event code ' num2str(evts(ind)) ' in trial #' ...
					num2str(t) ]);
		end
	end
end
events.success(isnan(events.success))=false;

