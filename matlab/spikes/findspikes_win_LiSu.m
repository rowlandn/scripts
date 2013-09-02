function spikes_SE = findspikes_win_LiSu(traces, fs, thres, varargin)
% 
% FINDSPIKES_WIN2_SE This function performs spike discrimination on  single
% peaks exceeding or within a threshold and/or time window range (ms). Peaks can
% be either positive or negative (defined by user) and do not have to have
% a  baseline of zero. Enter +1 for positive-going spikes and -1 for
% negative-going spikes.  Traces must be a matrix  with the rows
% representing voltages at successive sampling times  and the columns
% representing different traces.  The sampling frequency (Fs) in  kHz is
% used to scale the trace and spike times into ms. Output consists of a
% cell of spike time arrays with the row of the cell array representing the
% corresponding trace from which the spikes were discriminated. The second
% column of the cell array represents correspoding spike peak values.
%
%  spikes_SE = findspikes_win2_SE(traces, Fs, single_threshold)
%                       
%                            --- or ---
%
%  spikes_SE = findspikes_win2_SE(traces, Fs, [min_amp max_amp], direction, [min_time max_time], 'plot')
%
% 
%  Example 1: A = load_NS_SE('/Raw/bgv05/bgv0512a.data',1:10,2,10);
%             spikeTime = findspikes_win2_SE(A.Ch_2,10,-.3,1,'plot');
%             (to use a single amplitude threshold, use this notation)
%
%  Example 2: A = load_NS_SE('/Raw/bgv05/bgv0512a.data',1:10,2,10);
%             spikeTime = findspikes_win2_SE(A.Ch_2,10,[-.3 0],1,'plot');
%             (to use an amplitude threshold range, use this notation)
% 
%  Example 3: A = load_NS_SE('/Raw/bgv05/bgv0512a.data',1:10,2,10);
%             spikeTime = findspikes_win2_SE(A.Ch_2,10,[-.3 0],1,[.1 2],'plot');
%             (to use an amplitude threshold and time range in ms, use this notation)







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
            
%     if ~exist('plotit')
%         plotit = '';
%     end

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

% Set up waitbar
progbar = waitbar(0, 'Finding Spikes ...');

    % start to find spikes ===========================
    thres = thres*direction;
    thresh_min=min(thres);
    
    for idx = 1:size(traces,2)
        
        waitbar(idx/size(traces,2), progbar)

        trace = traces(:,idx);
        
        % flip the trace and threshold up-side-down to find the down-pointing spikes.
        trace = trace*direction;
        
        left_edges = find(trace(1:end-1) < thresh_min & trace(2:end) >= thresh_min) + 1;  % find rising slopes across the threshold
        right_edges = find(trace(1:end-1) >= thresh_min & trace(2:end) < thresh_min);  % find falling slopes across the threshold
        
        if right_edges(1) < left_edges(1)    % match the left and right edges of each window.
            right_edges(1)=[];
        end
        spike_num=min(length(left_edges), length(right_edges));
        
        % eliminate the time windows out of the range.
        if exist('win_range')
            data_range=win_range.*fs;
            win_width=right_edges(1:spike_num)-left_edges(1:spike_num);
            out_of_data_range_idx=find(win_width < data_range(1) | win_width > data_range(2));
            left_edges(out_of_data_range_idx)=[];
            right_edges(out_of_data_range_idx)=[];                
            spike_num=min(length(left_edges), length(right_edges));
        end
        
        timeidx=zeros(spike_num,1); peaks=zeros(spike_num,1);
        if spike_num > 0
            for sp=1:spike_num
                [peaks(sp) timeidx(sp)]=max(trace(left_edges(sp):right_edges(sp)));
            end
            timeidx=timeidx+left_edges(1:spike_num)-1;
	
            % eliminate the spikes that exceed the threshold range.
            if length(thres)>1 
                bigSpikes=find(peaks>max(thres));
                peaks(bigSpikes)=[];
                timeidx(bigSpikes)=[];
            end
        end
            assignin('base','timeidx',timeidx)
        spikes_SE{idx,1} = timeidx'/fs;
        peaks=peaks*direction;
        spikes_SE{idx,2} = peaks';
    end

    %Close Progess Bar
    close(progbar)

    % flip the traces back.
    thres=thres*direction;
    trace=trace*direction;
    
    if isequal(plotit,'plot')
        figure
        m_time = [1:size(trace,1)]'/(fs);
        plot(m_time,trace,'k'); hold on
        for I=1:length(thres)
            plot([0 m_time(end)],[thres(I) thres(I)],'b');
        end
        plot(timeidx/fs, peaks, 'ro'); hold off
        ylabel('V/uV')
        xlabel('t/ms')
        zoom on
    else
    end
return