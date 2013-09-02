function gauss = glr_SE(spike_times, gauss_width, Fs, max_time, plotit)

% glr_SE This function implements an algorithm known as  gaussian local
% rate coding (Paulin,M, 1996).  Individual  spikes in a spike train are
% convolved with a gaussian  curve, and the curves are then summed for each
% trace.   The resulting gaussian trace represents the instantaneous 
% frequency of the spike train in analog form, which can  then be
% correlated with other analog waveforms, such as  simultaneous recordings
% of EEG activity.  Spike times  should consist of a cell with individual
% rows representing  the spike times of individual traces (e.g., the output
% of the  findspikes_win_SE function.)  The user-defined width (ms) of 
% each gaussian curve is equal to the time that elapses  between the peak
% of the curve and the value 1/e^1/2.  Each  curve has a normalized area of
% 1. The sampling frequency  (Fs in Hz) is used to scale the time of the trace
% into a user-defined  maximum time, which is the desired length of the
% gaussian  trace(s) in ms (in order to perform auto- and
% cross-correlations,  both analog waveforms must be the same length.
% Therefore, the  user should match the length of the gaussian trace(s) to
% the  length of the other analog waveform(s).) Output of the function 
% consists of both individual traces and the average of the  individual
% traces. If plotit is 1, the average trace is plotted.
%
% gauss_trace = glr_SE(spike_times, gauss_width, Fs, max_time, plotit)
%  
% Example 1:  A = loadtraces_SE('/Raw/viv05/viv0515b.all','1-10',4); 
%             spike_times = findspikes_win_SE(A,10,{-100 -20 .1 1},1); 
%             gauss = glr_SE(spike_times,5,10,10000,1);


% % Scale parameters into data points
% gauss_width = 20;
% Fs = 10;
% max_time = 10000;
% plotit = 1;
% spike_times = a;

% Set up waitbar
progbar = waitbar(0, 'Convolving...');


for i = 1:size(spike_times,1)
    
    waitbar(i/length(spike_times), progbar)
    
    % Check desired length of signal
    if max(spike_times{i,1}) > max_time
    warndlg(['Desired maximum signal length too short.  Please',...
        ' enter a desired maximum signal length that exceeds the',...
        ' maximum spike time.']);
    break
    end
    
    % Scale parameters into data points
    spike_index = round(spike_times{i,1}*Fs);
    max_index = max_time*Fs;
    spike_train = zeros(max_index,1);
    spike_train(spike_index) = 1;


    % Calculate and scale gaussian curve
    sigma = gauss_width*Fs*10;
    time_range = [-100*gauss_width:100*gauss_width]*10;
    %gauss_curve = 1./(sqrt(2*pi*sigma))*exp(-time_range.^2/(2*sigma^2));
    gauss_curve = (1./(sqrt(2*pi)*sigma))*(exp(-time_range.^2/(2*sigma^2)));
    gauss_curve = gauss_curve/trapz(gauss_curve);
    
    glr = conv(spike_train, gauss_curve)*10000;
    m_cutoff = fix((length(glr)-length(spike_train))/2);
    glr = glr(1+m_cutoff:max_index+m_cutoff);
    
    gauss_mean(:,i) = glr;
    
    gauss.ind_traces(:,i) = glr;
    
    
end

% Scale time axis of individual traces
gauss.ind_traces(:,size(gauss.ind_traces,2)+1) = [1:size(gauss.ind_traces,1)]'/Fs;

% Calculate average of individual traces
mean_gauss = mean(gauss_mean,2);
mean_gauss(:,2) = [1:size(mean_gauss,1)]'/Fs;
gauss.avg_trace = mean_gauss;

% Close Progess Bar
close(progbar)

if plotit==1
    figure;
    plot(gauss.avg_trace(:,2),gauss.avg_trace(:,1));
    xlabel('ms')
    max_gauss_trace = max(gauss.avg_trace(:,1));
    axis([0 max_time 0 max_gauss_trace+.2*max_gauss_trace])
end

