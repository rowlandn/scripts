function emg_offset = DetectEMG_off(time, emg, epoch)
%   Detects onset of emg, created for use in time_psd analysis
%	Input:
%       time    - vector containing time, assuming sampling rate of 1kHz
%       emg     - vector containing emg data, assuming sampling rate of
%       1kHz
%       epoch   - rest epoch timestamps to be used for emg alignment
%	Output:
%       emg_offset   - vector containing user-defined emg offset timestamps
%                   
%   Adapted from DetectEMG by: A.Crowell 12/14/2009

% Initialize variables
Fs = 1000;  % emg sampling rate
THRESH = 7;	% percent of cntl SD
EMG_bl = [2 3];    % seconds before epoch timestamp used for baseline EMG activity
EMG_test = [-1 1]; % seconds around epoch timestamp to test EMG
n_epoch = length(epoch);
emg_offset = zeros(1,n_epoch);
c = [];

h = figure;
for i = 1:n_epoch
	cntl_inds = find(time>(epoch(i)+EMG_bl(1)) & time<(epoch(i)+EMG_bl(2)));
	test_inds = find(time>(epoch(i)+EMG_test(1)) & time<(epoch(i)+EMG_test(2)));
	cntl_mn = mean(emg(cntl_inds));
	cntl_sd = mad(emg(cntl_inds)); %mad = mean/median abs dev
	
	start_ind = min(find(abs(emg(test_inds))> THRESH*cntl_sd )); %#ok<MXFND>
	start_time = time(test_inds(start_ind));
	
	plot(time,emg,'k')
	xlim([epoch(i)-5 epoch(i)+5]);
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
    ylabel('EMG activity');
	
    if i==1
        disp(['If offset is correct, hit return to move to next epoch.' sprintf('\n')...
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
	emg_offset(i) = start_time;
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







    

