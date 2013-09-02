function autocorr = auto_corr_SE(X,Fs,maxlag,sp1,sp2,sp3)

% auto_corr_SE This function performs an auto-correlation of a single signal.  
% Individual traces within the signal should be represented as column 
% vectors. Sampling frequency (Fs) should be given in kHz and is used 
% along with maxlag (s) to scale the trace into seconds. If multiple 
% auto-correlograms are to be placed on the same page, subplot (sp) 
% coordinates can be given.  If not, these can be omitted.  
%
% autocorr = auto_corr_SE(Signal,Fs,maxlag,sp1,sp2,sp3)
%
% Example 1: A = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',1);
%            autocorr = auto_corr_SE(A,10,2);
%
% Example 2: A = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',1);
%            autocorr = auto_corr_SE(A,10,2,2,2,1);
%                 (if you want to specify a certain
%                  sublot, use this notation)
%
% Example 3: Neuron = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',4);
%            spike_times = findspikes_win_SE(Neuron,10,-200,-40,.1,1,0);
%            gauss = glr_SE(spike_times,5,10,10000,0);
%            autocorr = auto_corr_SE(gauss.ind_traces(:,1:10),10,2);
%                    (to correlate spike times transformed with the glr_SE function
%                    and another analog signal, use this notation)

Signal_1 = X;
Signal_2 = X;

maxlag = maxlag*Fs*1000;

% Concatenate individual traces (represented as columns) in matrix into a
% single array
cat_Signal_1 = reshape(Signal_1,prod(size(Signal_1)),1);
cat_Signal_2 = reshape(Signal_2,prod(size(Signal_2)),1);
[cov,lags] = xcov(cat_Signal_1,cat_Signal_2,maxlag,'coeff');

cov_null = xcov(cat_Signal_1,flipud(cat_Signal_2),maxlag,'coeff');
CI_cov = 3.1*std(cov_null);

if nargin == 6
    subplot(sp1,sp2,sp3)
else
end 
CI_bar = patch(maxlag/(Fs*1000)*[-1 -1 1 1],CI_cov*[-1 1 1 -1],'y');
set(CI_bar,'facecolor',[.85 .85 .85],'edgecolor',[1 1 1])
hold on
cov = cov.';
% assignin('base','cov',cov)
% assignin('base','lags',lags)
plot(lags/(Fs*1000),cov,'k','LineWidth',2); hold on
ylabel('Correlation');
xlabel('Time Lag (s)');
y_lim = get(gca,'ylim');
x_lim = get(gca,'xlim');

axis([x_lim(1)+.2*x_lim(1) x_lim(2)+.2*x_lim(2) -1.2 1.2])

plot([0 0],[-1.2 1.2],'k')

    

[max_cov,index_max_cov] = max(cov);
[min_cov,index_min_cov] = min(cov);
% assignin('base','max_cov',max_cov)
% 
% assignin('base','min_cov',min_cov)
% assignin('base','index_max_cov',index_max_cov)
% 
% assignin('base','index_min_cov',index_min_cov)

if abs(max_cov) > abs(min_cov) & round((index_max_cov-maxlag-1))/(Fs*1000) >= 0
    text(5000,.4,num2str(round((index_max_cov-maxlag-1))/(Fs*1000)))
    peak_value = max_cov;
    peak_index = round((index_max_cov-maxlag-1))/(Fs*1000);
elseif abs(max_cov) > abs(min_cov) & round((index_max_cov-maxlag-1))/(Fs*1000) < 0
    text(-5000,.4,num2str(round((index_max_cov-maxlag-1))/(Fs*1000)))
    peak_value = max_cov;
    peak_index = round((index_max_cov-maxlag-1))/(Fs*1000);
elseif abs(min_cov) > abs(max_cov) & round((index_min_cov-maxlag-1))/(Fs*1000) >= 0
    text(5000,-.4,num2str(round((index_min_cov-maxlag-1))/(Fs*1000)))
    peak_value = min_cov;
    peak_index = round((index_min_cov-maxlag-1))/(Fs*1000);
elseif abs(min_cov) > abs(max_cov) & round((index_min_cov-maxlag-1))/(Fs*1000) < 0
    text(-5000,-.4,num2str(round((index_min_cov-maxlag-1))/(Fs*1000)))
    peak_value = min_cov;
    peak_index = round((index_min_cov-maxlag-1))/(Fs*1000);
end

% assignin('base','peak_value',peak_value)
% 
% assignin('base','peak_index',peak_index)
peak_latency = peak_index;

autocorr.cov = cov;
autocorr.lags = lags/(Fs*1000);
autocorr.CI_bar_x = get(CI_bar,'XData');
autocorr.CI_bar_y = get(CI_bar,'YData');
autocorr.peak_value = peak_value;
autocorr.peak_latency = peak_index;

