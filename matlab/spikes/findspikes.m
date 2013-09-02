function [spikeTime, spikePeak] = findspikes(traces, fs, thres, varargin)
% 
% FINDSPIKES2_SE This function performs spike discrimination on 
% single peaks exceeding or within a threshold and/or time window range.
%
% Syntax:
%
%  spikeTime = findspikes2_SE(traces, fs, threshold)
%  [spikeTime spikePeak] = findspikes2_SE(traces, fs, threshold [,direction] [,win_range] [,'plot'])
% 
% Description:
%
%  traces   : Multiple traces of signal. each trace in a column
%  fs       : Sampling frequency, in KHz
%  threshold: Either a scalar, or [thres1 thres2] to define a range
%  direction: Optional. 
%             A positive number to find positive-going spikes, and vice versa. 
%             The default value is +1 when threshold is a scalar,
%             is sign(thres2-thres1) when threshold is a vector.
%  win_range: Optional. 
%             [win_min win_max] to define a time window range. in ms.
%  'plot'   : Optional. Plot the result. Not to plot by default.
%  spikeTime: Returns a cell array, each cell is a vector of spike times in each trial.
%  spikePeak: Returns the peak values of each spikes.
% 
% Samples:
%
%  spikeTime=findspikes2_SE(signal, 10, 0.25);
%  spikeTime=findspikes2_SE(signal, 10, 0.25, 'plot');
%  spikeTime=findspikes2_SE(signal, 10, [-0.2 -0.5]);
%  spikeTime=findspikes2_SE(signal, 10, [-0.5 -0.2], -1);
%  spikeTime=findspikes2_SE(signal, 10, [-0.2 -0.3], [0.1 2], 'plot');

    % assign the arguments========================
    error(nargchk(3,6,nargin))
    
    for k=1:nargin-3
        if ischar(varargin{k})
            plotit=varargin{k};
        elseif isnumeric(varargin{k})
            if length(varargin{k})==1
                direction=varargin{k};
            else
                win_range=varargin{k};
            end
        else
            error('error')
        end
    end
            
    if ~exist('plotit')
        plotit = '';
    end

    if ~exist('direction')
        if length(thres)==1
            direction = 1;
        else
            direction=thres(2)-thres(1);
        end
    end

    direction=sign(direction); 
    if direction==0
        error('threshold range or direction is zero')
    end

    % start to find spikes ===========================
    thres = thres*direction;
    thresh_min=min(thres);
    
    for idx = 1:size(traces,2)
        
        trace = traces(:,idx);
        
        % flip the trace and threshold up-side-down to find the down-pointing spikes.
        trace = trace*direction;
        
        left_edges = find(trace(1:end-1) < thresh_min & trace(2:end) >= thresh_min);  % find rising slopes across the threshold
        right_edges = find(trace(1:end-1) >= thresh_min & trace(2:end) < thresh_min);  % find falling slopes across the threshold
        
		if ~(isempty(left_edges) | isempty(right_edges))
            if right_edges(1) < left_edges(1)    % match the left and right edges of each window.
                right_edges(1)=[];
            end
		end
        spike_num=min(length(left_edges), length(right_edges));
		left_edges=left_edges(1:spike_num);
		right_edges=right_edges(1:spike_num);
		
        % eliminate the time windows out of the range.
        if spike_num > 0 & exist('win_range')
            data_range=win_range.*fs;
			left_time=left_edges + (thresh_min - trace(left_edges)) ./ (trace(left_edges+1) - trace(left_edges));
			right_time=right_edges + (thresh_min - trace(right_edges)) ./ (trace(right_edges+1) - trace(right_edges));
            win_width=right_time-left_time;
            out_of_data_range_idx=find(win_width < data_range(1) | win_width > data_range(2));
            left_edges(out_of_data_range_idx)=[];
            right_edges(out_of_data_range_idx)=[];                
            spike_num=min(length(left_edges), length(right_edges));
        end

		% merge the neighbor windows too close to each other.
		refractory=1;
		if spike_num > 0 
			interval=left_edges(2:end)-right_edges(1:end-1);
			smallInterval=find(interval < refractory*fs);
			right_edges(smallInterval)=[];
			left_edges(smallInterval+1)=[];
            spike_num=min(length(left_edges), length(right_edges));
		end
		
        timeidx=zeros(spike_num,1); peaks=zeros(spike_num,1);
        if spike_num > 0
            for sp=1:spike_num
                [peaks(sp) timeidx(sp)]=max(trace(left_edges(sp):right_edges(sp)));
            end
            timeidx=timeidx+left_edges-1;
	
            % eliminate the spikes exceed the max threshold.
            if length(thres)>1 
                bigSpikes=find(peaks>max(thres));
                peaks(bigSpikes)=[];
                timeidx(bigSpikes)=[];
            end
			
% 			% eliminate small intervals.
% 			interval=diff(timeidx)/fs;
% 			smallInterval=find(interval<1);
% 			rmlist=[];
% 			for k=1:length(smallInterval)
% 				if min(trace(timeidx(smallInterval(k)):timeidx(smallInterval(k)+1)))>thresh_min/2
% 					rmlist(end+1)=smallInterval(k);
% 					if peaks(rmlist(end))>peaks(rmlist(end)+1)
% 						rmlist(end)=rmlist(end)+1;
% 					end
% 				end
% 			end
% 			peaks(rmlist)=[];
% 			timeidx(rmlist)=[];
			
			spike_num=length(timeidx);
        end
            
        spikeTime{idx} = timeidx/fs;
        if nargout==2
            spikePeak{idx} = peaks;
        end
    end
    spikeTime = spikeTime';
    spikePeak = spikePeak';
    
    % flip the traces back.
    peaks=peaks*direction;
    thres=thres*direction;
    trace=trace*direction;
    
    if isequal(plotit, 'plot')
        figure
        m_time = [1:size(trace,1)]'/(fs);
        plot(m_time,trace,'k'); hold on
        for I=1:length(thres)
            plot([0 m_time(end)],[thres(I) thres(I)],'b');
        end
        plot(timeidx/fs, peaks, 'ro'); 
		hold off
        ylabel('V/uV')
        xlabel('t/ms')
        zoom on
    end
return