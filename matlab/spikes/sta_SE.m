function sta = sta_SE(Analog_Signal,spike_times,lag,Fs,plotit,sp1,sp2,sp3)

% sta_SE This function computes the spike-triggered average of an analog
% signal around given spike times.  Lag (ms) defines the time of the analog
% signal before and after each spike time. Spike times should be supplied
% as a cell containing rows of spike time arrays which represent individual
% traces (e.g., the SPIKE_TIMES variable from the findspikes_win_SE
% function). Sample rate should  be given in kHz.   The spike times are
% also randomly shuffled 50 times to produce 99% confidence intervals for
% the original spike-triggered average.  If  multiple spike-triggered
% averages are to be placed on the same page, subplot (sp)  coordinates can
% be given.  If not, these can be omitted.
% 
%  sta = sta_SE(Analog_Signal,spike_times,lag,Fs,sp1,sp2,sp3)
% 
%  Example 1: Analog_Signal = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',1);
%             Neuron = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',4);
%             spike_times = findspikes_win_SE(Neuron,10,-200,-40,.1,1,0);
%             sta = sta_SE(Analog_Signal,spike_times,500,10,1);
%              
%  Example 2: Analog_Signal = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',1);
%             Neuron = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',4);
%             spike_times = findspikes_win_SE(Neuron,10,-200,-40,.1,1,0);
%             sta = sta_SE(Analog_Signal,spike_times,500,10,1,2,2,1);
%                               (if you want to specify a certain
%                                sublot, use this notation)

              
% % % % Debugging code; comment out
% Analog_Signal = loadtraces_SE('/F/Data/Raw/viv05/viv0503b.all','1-10',1);
% DCN = loadtraces_SE('/F/Data/Raw/viv05/viv0503b.all','1-10',4);
% spike_times = findspikes_win_SE(DCN,10,-200,-40,.1,1,0);
% Fs = 10;
% lag = 500;
% plotit = 1;

progbar = waitbar(0, 'Computing average with original spike times...');
%%%%% Here we find the spike-triggered average for the %%%%%%%%%%%%
%%%%% original signal and spike times.                 %%%%%%%%%%%%
 

for i = 1:size(spike_times,1)
    
    waitbar(i/size(spike_times,1), progbar)
    spike_index = round(spike_times{i,1}*Fs);
    spike_train = zeros(size(Analog_Signal(:,i)));
    spike_train(spike_index) = 1;

    [orig_sta(:,i) orig_m_lags] = xcorr(Analog_Signal(:,i),spike_train(:,1),lag*Fs);
    
    orig_sta(:,i) = orig_sta(:,i)/length(spike_times{i,1});
    
end
close(progbar)

mean_orig_sta = mean(orig_sta,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear spike_index spike_train


%%%%% Here, we find new, shuffled spike times for each trace %%%%% 
%%%%% that we will use to make new spike-triggered averages. %%%%% 
%%%%% The new, shuffled spike-triggered averages will be used %%%%
%%%%% to calculate confidence intervals for the original spike- %%
%%%%% triggered average.                                     %%%%%
progbar = waitbar(0, 'Now, shuffling spike times...');
for m = 1:50
        
        %assignin('base','m',m)

        waitbar(m/size(spike_times,1), progbar)
        for l = 1:size(spike_times,1)
            isi_str = ['isi_trace_',num2str(l),' = diff(spike_times{',num2str(l),',1});'];
            eval(isi_str)
            randperm_str_1 = ['randperm_size_trace_',num2str(l),' = randperm(size(spike_times{',num2str(l),',1},2)-1);'];
            eval(randperm_str_1)
            shuf_str_1 = ['shuf_isis_trace_',num2str(l),' = isi_trace_',num2str(l),'(randperm_size_trace_',num2str(l),');'];
            eval(shuf_str_1);       
            shuf_str_2 = ['shuf_spike_times_',num2str(m),'{',num2str(l),',1} = [1 (cumsum(shuf_isis_trace_',num2str(l),')+1)];'];            
            eval(shuf_str_2)
        end
        
        spike_times_shuf_str_2 = ['spike_times_shuf = shuf_spike_times_',num2str(m),';'];
        eval(spike_times_shuf_str_2)
        
end
close(progbar)

progbar = waitbar(0, 'Now, computing average with shuffled spike times...');
for i = 1:size(spike_times_shuf,1)

    waitbar(i/size(spike_times,1), progbar)
    spike_index = round(spike_times_shuf{i,1}*Fs);
    spike_train = zeros(size(Analog_Signal(:,i)));
    spike_train(spike_index) = 1;

    [shuf_sta(:,i) shuf_m_lags] = xcorr(Analog_Signal(:,i),spike_train(:,1),lag*Fs);
    
    shuf_sta(:,i) = shuf_sta(:,i)/length(spike_times_shuf{i,1});
    
end
close(progbar)

mean_shuf_sta_y = mean(shuf_sta,2);
mean_shuf_sta_x = shuf_m_lags'/Fs;
std3_plus_mean_shuf_sta_y = mean_shuf_sta_y + 3.1*std(mean_shuf_sta_y);
std3_minus_mean_shuf_sta_y = mean_shuf_sta_y - 3.1*std(mean_shuf_sta_y);
mean_std3_plus_mean_shuf_sta_y = mean(std3_plus_mean_shuf_sta_y);
mean_std3_minus_mean_shuf_sta_y = mean(std3_minus_mean_shuf_sta_y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Plot  

%if plotit == 1 
    if nargin == 8
        subplot(sp1,sp2,sp3)
    else
    end
    m_t = orig_m_lags/Fs;
    CI_bar_offset = patch(lag*[-1 -1 1 1],...
             [mean_std3_plus_mean_shuf_sta_y mean_std3_minus_mean_shuf_sta_y mean_std3_minus_mean_shuf_sta_y mean_std3_plus_mean_shuf_sta_y],...
             'w'); hold on
    set(CI_bar_offset,'facecolor',[.85 .85 .85],'edgecolor',[1 1 1])
     plot(m_t, mean_orig_sta,'k'); hold on
    %plot(mean_shuf_sta_x,mean_shuf_sta_y)
%     plot(mean_shuf_sta_x,std3_plus_mean_shuf_sta_y,'r')
%     plot(mean_shuf_sta_x,std3_minus_mean_shuf_sta_y,'r')    
    x_lim = get(gca,'xlim');
    y_lim = get(gca,'ylim');
    plot([0 0],[y_lim(1) y_lim(2)],'k')
    
    ylabel('Spike Triggered Average')
    xlabel('Lag (ms)')
    %end


sta.sta_y = mean_orig_sta;
sta.sta_x = orig_m_lags'/Fs;
sta.std_plus_shuf_mean = std3_plus_mean_shuf_sta_y;
sta.std_minus_shuf_mean = std3_minus_mean_shuf_sta_y;
sta.CI_bar_x = get(CI_bar_offset,'XData');
sta.CI_bar_y = get(CI_bar_offset,'YData');

[min_sta_y,min_sta_i] = min(sta.sta_y);
[max_sta_y,max_sta_i] = max(sta.sta_y);

if abs(min_sta_y) > abs(max_sta_y)
    sta.sta_peak_lag = sta.sta_x(min_sta_i);
    sta.sta_peak_sta = min_sta_y;
elseif abs(max_sta_y) > abs(min_sta_y)
    sta.sta_peak_lag = sta.sta_x(max_sta_i);
    sta.sta_peak_sta = max_sta_y;
end

hold on; plot(sta.sta_peak_lag,sta.sta_peak_sta,'ro')

