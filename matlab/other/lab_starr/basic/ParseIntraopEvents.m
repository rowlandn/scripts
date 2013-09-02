function	events = ParseIntraopEvents(time, evt_chan)
% function	events = ParseIntraopEvents(time, evt_chan)
%
%	Detects event times in event channels
%	Assumes 1 kHz sampling rate for 
%	Input:
%		time - vector of time
%		evt_chan - analog channel in which voltages reflect 
%					switch & light combinations
%
%	Output:
%		events - structure containing vectors of event times 
%				
%	RST - 2006-04-05

MARGIN = 0.05; 
OFF_NONE =  [5.14-MARGIN 5.14+MARGIN];  
HOME_NONE = [4.37-MARGIN 4.37+MARGIN];
HOME_LEFT = [3.15-MARGIN 3.15+MARGIN];
OFF_LEFT = [3.54-MARGIN 3.54+MARGIN];
HOME_RIGHT = [1.57-MARGIN 1.57+MARGIN];
OFF_RIGHT = [1.67-MARGIN 1.67+MARGIN];
%______________________________
% Commented out by Sho Shimamoto on 6/27/07 for task hardware troubleshoot
% Task voltage levels were decreased by a factor of 5
% MARGIN = 0.005;
% OFF_NONE =  [1.028-MARGIN 1.028+MARGIN];
% HOME_NONE = [0.874-MARGIN 0.874+MARGIN];
% HOME_LEFT = [0.63-MARGIN 0.63+MARGIN];
% OFF_LEFT = [0.707-MARGIN 0.707+MARGIN];
% HOME_RIGHT = [0.313-MARGIN 0.313+MARGIN];
% OFF_RIGHT = [0.330-MARGIN 0.330+MARGIN];
%______________________________


[home_ev,home_val] = getevts(time, evt_chan, NaN, HOME_NONE);

[left_ev,left_val] = getevts(time, evt_chan, HOME_NONE, HOME_LEFT);
[left_mv_ev,left_mv_val] = getevts(time, evt_chan, HOME_LEFT, OFF_LEFT);
[right_ev,right_val] = getevts(time, evt_chan, HOME_NONE, HOME_RIGHT);
[right_mv_ev,right_mv_val] = getevts(time, evt_chan, HOME_RIGHT, OFF_RIGHT);

h = figure;
plot(time,evt_chan,'k')
axis tight
ylm = ylim;
hold on
plotevts( home_ev, home_val, ylm(2), 'c', 'Hm');

plotevts( left_ev, left_val, ylm(2)-0.15, 'r', 'LftLt');
plotevts( left_mv_ev, left_mv_val, ylm(2)-.3, 'm', 'LfMv');
plotevts( right_ev, right_val, ylm(2)-0.15, 'g', 'RtLt');
plotevts( right_mv_ev, right_mv_val, ylm(2)-.3, 'b', 'RtMv');

display('Please review event parsing.  Hit any key to continue.');
pause

% Assemble into coherent structure & check for errors
left_events = CompileEvents(home_ev,left_ev,left_mv_ev,...
	[right_ev right_mv_ev],'L');
right_events = CompileEvents(home_ev,right_ev,right_mv_ev,...
	[left_ev left_mv_ev],'R');

% Interleave Right & Left trials into one structure
events.home = [left_events.home right_events.home];
events.light = [left_events.light right_events.light];
events.move = [left_events.move right_events.move];
events.target = [left_events.target right_events.target];

[events.light,inds] = sort(events.light);
events.home = events.home(inds);
events.move = events.move(inds);
events.target = events.target(inds);
