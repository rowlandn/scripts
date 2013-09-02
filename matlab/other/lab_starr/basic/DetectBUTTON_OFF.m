function signal_onset = DetectBUTTON_OFF(time, signal,signal2, epoch)
%   Detects onset of signal, created for use in time_psd analysis
%	Input:
%       time    - vector containing time, assuming sampling rate of 1kHz
%       signal     - vector containing signal data, assuming sampling rate of 1kHz
%       epoch   - active epoch timestamps to be used for signal alignment
%	Output:
%       signal_onset   - vector containing user-defined signal onset timestamps
%                   
%   Created by: S.Shimamoto 11/6/2008

% Initialize variables
THRESH = min(signal(1000:end-1000))+(max(signal(1000:end-1000))-min(signal(1000:end-1000)))/2;%7;	% percent of cntl SD
signal_bl = [0 3];    % seconds before epoch timestamp used for baseline signal activity
signal_test = [0 5]; % seconds around epoch timestamp to test signal
n_epoch = length(epoch);
signal_onset = zeros(1,n_epoch);
c = [];

h = figure;
for i = 1:n_epoch
	cntl_inds = find(time>(epoch(i)+signal_bl(1)) & time<(epoch(i)+signal_bl(2)));
	test_inds = find(time>(epoch(i)+signal_test(1)) & time<(epoch(i)+signal_test(2)));
    cntl_mn = mean(signal(cntl_inds));
	cntl_sd = mad(signal(cntl_inds));
	
	start_ind = min(find(abs(signal(test_inds))< THRESH )); %#ok<MXFND>
	start_time = time(test_inds(start_ind));
	
	plot(time,signal,'k')
    hold on
    plot(time,signal2,'b')
	xlim([epoch(i)-10 epoch(i)+10]);
	xlm = xlim;
	hold on
	plot(xlm,...
		[cntl_mn-THRESH*cntl_sd cntl_mn-THRESH*cntl_sd],'r:');
	plot(xlm,...
		[cntl_mn+THRESH*cntl_sd cntl_mn+THRESH*cntl_sd],'r:');
    ylm = ylim;
    plot([epoch(i) epoch(i)],ylm,'k:');

	if ~isempty(start_time)
		c = plotevts( start_time, ylm(1), ylm(2), 'r');
    end

    title(['epoch #' num2str(i)]);
    xlabel('time (sec)');
    ylabel('signal activity');
	
    if i==1
        disp(['If onset is correct, hit return to move to next epoch.' sprintf('\n')...
            'Otherwise, correct onset as needed or click outside box to delete, then hit return.']);
    end
	[x,y] = ginput(1);
	while ~isempty(x)
		if ishandle(c)
			delete(c)
		end
		if inrange(x,xlm) && inrange(y,ylm)
			start_time = x;
			c = plotevts( start_time, ylm(1), ylm(2), 'r');
		else
			start_time = [];
		end
		[x,y] = ginput(1);
	end
	hold off
	if isempty(start_time)
		start_time=NaN;
	end
	move_onset(i) = start_time;
end

close(h)

return;

% subfunction
function in = inrange(chan, range)
% function in = inrange(chan, range)
%
%	Return logic if chan is within range
%
%	RST 2006-04-05

in = chan > range(1) & chan < range(2);







    

