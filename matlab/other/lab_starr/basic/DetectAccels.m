function	events = DetectAccels(time, accel, events, logic)
% function	events = DetectAccels(time, accel, events, logic)
%
%	Detects 
%	Assumes 1 kHz sampling rate for 
%	Input:
%		time - vector of time
%		accel - analog channel in which voltages reflect 
%					arm acceleration
%		events - structure of trials
%       logic - '1' indicates accel data, '0' indicates emg data
%	Output:
%		events - structure containing vectors of event times 
%				
%	RST - 2006-04-05
%   Edited 1-22-2008 by SS to detect both accel and emg data


THRESH1 = 7;	% percent of cntl SD
THRESH2 = 2;	% percent of cntl SD
n_trials = length(events.light);
c = [];

h = figure;
for i = 1:n_trials
	hm = events.home(i);
	lt = events.light(i);
	mv = events.move(i);
	if isnan(hm) || isnan(lt) || isnan(mv)
		events.accel(i) = NaN;
		continue
	end
	cntl_inds = find(time>hm & time<lt);
	test_inds = find(time>lt & time<mv);
	cntl_mn = mean(accel(cntl_inds));
	cntl_sd = mad(accel(cntl_inds));
	
	start_ind = min(find( abs(accel(test_inds))> THRESH1*cntl_sd ));
	start_time = time(test_inds(start_ind));
	
	plot(time,accel,'k')
	xlim([events.home(i) events.move(i)+1]);
	xlm = xlim;
	ylm = ylim;
	hold on
	plot([events.home(i) events.move(i)+1],...
		[cntl_mn-THRESH1*cntl_sd cntl_mn-THRESH1*cntl_sd],'r:');
	plot([events.home(i) events.move(i)+1],...
		[cntl_mn+THRESH1*cntl_sd cntl_mn+THRESH1*cntl_sd],'r:');
	
	plotevts( events.light, ylm(1), ylm(2)-.15, 'r', 'Light');
	plotevts( events.move, ylm(1), ylm(2)-.15, 'm', 'Switch');

	if ~isempty(start_time)
		c = plotevts( start_time, ylm(1), ylm(2), 'b');
	end

	disp('Correct onset if desired. Click outside box to delete. Otherwise hit return.');
	[x,y] = ginput(1);
	while ~isempty(x)
		if ishandle(c)
			delete(c)
		end
		if inrange(x,xlm) && inrange(y,ylm)
			start_time = x;
			c = plotevts( start_time, ylm(1), ylm(2), 'b');
		else
			start_time = [];
		end
		[x,y] = ginput(1);
	end
	hold off
	if isempty(start_time)
		start_time=NaN;
	end
	events.accel(i) = start_time;
end

close(h)

if logic  % check to see if data is accel 
    return;
else % if data is emg, store in new field called 'emg' then delete 'accel'
    events.emg = events.accel;
    events = rmfield(events,'accel');
end
return;

    

