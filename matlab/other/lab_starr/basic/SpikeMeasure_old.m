function spk = SpikeMeasure( spk_wf, fs, fname, unitname, show)
% spk = SpikeMeasure( spk_wf, fs, unitname, show )
% This program computes measures of spike snippets
% Inputs:
%	spk.wf = array of spike snippets (x samples long -by- y snippets)
%	fs = sampling frequency of snippets
%	fname - file being processed
%	unitname = name of snippets from NEX file
%	show - optional flag to control plotting (default=true)
%
% Outputs:
%	spk = structure containing results
%
% Created by RST, 2005-08-02
%
THRESH_PCT = 0.1;	%Percent max change used to detect onset/offset

error = false;

if ~exist('show','var')
	show = true;
end

dt = 1000/fs;	% Sample interval in msec

[dat_len,spk.n] = size(spk_wf);
spk.fname = fname;
spk.unit = unitname;
spk.mean = mean(spk_wf,2);

% Remove offset if present - seems to hurt more than help
%spk.mean = spk.mean - mean(spk.mean(1:3));	

% Find key points in mean spike
[spk.min, min_ind] = min(spk.mean);
[spk.max, max_ind] = max( spk.mean(min_ind:end) );
max_ind = max_ind + min_ind-1;

% Major negative phase almost always comes first
if max_ind > min_ind
	start_ind = min( find(spk.mean <= spk.min*THRESH_PCT) ) - 1;
	stop_ind = max( find(spk.mean >= spk.max*THRESH_PCT) );
	if stop_ind == dat_len
		spk.report = 'long';
	else
		spk.report = 'good';
		stop_ind = stop_ind+1;
	end
	
	if start_ind >1
	
		time = dt.*((1:dat_len) - start_ind);

		% Compute summary stats
		spk.start2max = time(max_ind);
		spk.start2stop = time(stop_ind);
		spk.min2max = time(max_ind) - time(min_ind);
		spk.area = sum( abs( spk.mean(start_ind:stop_ind) ));

	else
		error = true;
	end
else
	error = true;
end

if error
	time = dt.*(1:dat_len);
	spk.start2max = NaN;
	spk.start2stop = NaN;
	spk.min2max = NaN;
	spk.area = NaN;
	warning('Unable to measure this spike - wierd shape.');
	spk.report = 'error';
end

% Plot results
if show
	figure
	plot(time,spk.mean,'k','LineWidth',2);
	ysc = ylim;
	hold on
	axis tight
	xsc = xlim;
	ylim(ysc);
	if ~error
		plot([time(start_ind) time(start_ind)],ysc,'r');
		plot([time(stop_ind) time(stop_ind)],ysc,'r');
		plot(time(min_ind),spk.min,'ro','LineWidth',2);
		plot(time(max_ind),spk.max,'ro','LineWidth',2);
	else
		x = diff(xsc)/2 + xsc(1);
		y = diff(ysc)/2 + ysc(1);
		text( x, y, 'ERROR: Wierd Spike','FontSize',14,...
			'Color','r','HorizontalAlignment','center');
	end
	xlabel('msec');
	set(gca,'FontSize',8);
	title([fname '-:' unitname ],'Interpreter','none','FontSize',10);
end
