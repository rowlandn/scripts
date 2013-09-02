function crosscorr = cross_corr_SE(S,T,Fs,maxlag,threshold,sp1,sp2,sp3)

% function crosscorr = cross_corr_SE(S,T,Fs,maxlag,threshold,sp1,sp2,sp3)
%
% cross_corr_SE This function performs a cross-correlation between 2 
% signals - a source and a target signal. The target trace is slid across
% the source trace to calculate the correlation.  When there is a positive
% phase lag, it means the source leads the target and vice versa.  The
% signals  should be represented as equally sized  column vectors, with the
% rows representing voltages at successive sampling  times and different
% columns representing different traces. Simultaneously  acquired signals
% should have the same column index, such that S(:,i) will be  correlated
% with T(:,i).  Sampling frequency (Fs) should be given in kHz and is used
% along with maxlag (s) to scale the trace into seconds. Enter threshold
% for local peak discrimination. If  multiple cross-correlograms are to be
% placed on the same page, subplot (sp)   coordinates can be given.  If
% not, these can be omitted.  
%
% crosscorr = cross_corr_SE(S,T,Fs,maxlag,threshold,sp1,sp2,sp3)
%
% Example 1: A = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',1);
%            B = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',2);
%            crosscorr = cross_corr_SE(A,B,10,2);
%
% Example 2: A = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',1);
%            B = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',2);
%            crosscorr = cross_corr_SE(A,B,10,2,2,2,1);
%                 (if you want to specify only certain columns of a matrix
%                  to be correlated and you want to specify a certain
%                  sublot, use this notation)
%
% Example 3: Analog_Signal = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',1);
%            Neuron = load_PCDX_SE('/Raw/viv05/viv0503b.all','1-10',4);
%            spike_times = findspikes_win_SE(Neuron,10,-200,-40,.1,1,0);
%            gauss = glr_SE(spike_times,5,10,10000,0);
%            crosscorr = cross_corr_SE(Analog_Signal,gauss.ind_traces(:,1:10),10,2);
%                    (to correlate spike times transformed with the glr_SE function
%                    and another analog signal, use this notation)

Signal_1 = T;
Signal_2 = S;

maxlag = maxlag*Fs*1000;

% Concatenate individual traces (represented as columns) in matrix into a
% single array
cat_Signal_1 = reshape(Signal_1,prod(size(Signal_1)),1);
cat_Signal_2 = reshape(Signal_2,prod(size(Signal_2)),1);
[cov,lags] = xcov(cat_Signal_1,cat_Signal_2,maxlag,'coeff');
%assignin('base','lags',lags)
cov_null = xcov(cat_Signal_1,flipud(cat_Signal_2),maxlag,'coeff');
CI_cov = 3.1*std(cov_null);

if nargin == 8
    subplot(sp1,sp2,sp3)
else
end

CI_bar = patch(maxlag/(Fs*1000)*[-1 -1 1 1],CI_cov*[-1 1 1 -1],'y');
set(CI_bar,'facecolor',[.85 .85 .85],'edgecolor',[1 1 1])
hold on
cov = cov.';
plot(lags/(Fs*1000),cov,'k','LineWidth',2); hold on
%plot(lags/(Fs*1000),cov,'k','LineWidth',2); hold on
y_lim = get(gca,'ylim');
x_lim = get(gca,'xlim');

axis([x_lim(1)+.2*x_lim(1) x_lim(2)+.2*x_lim(2) y_lim(1)+.2*y_lim(1) y_lim(2)+.2*y_lim(2)])

%plot([0 0],[y_lim(1)+.2*y_lim(1) y_lim(2)+.2*y_lim(2)],'k')

plot([0 0],[-1 1],'k')
ylabel('Correlation');
xlabel('Time Lag (s)');


 
        


% [max_cov,index_max_cov] = max(cov);
% [min_cov,index_min_cov] = min(cov);
% 
% % assignin('base','max_cov',max_cov)
% % assignin('base','min_cov',min_cov)
% % assignin('base','index_max_cov',index_max_cov)
% % assignin('base','index_min_cov',index_min_cov)
% 
% 
% if abs(max_cov) > abs(min_cov) & round((index_max_cov-maxlag-1))/(Fs*1000) >= 0
%     text(5000,.4,num2str(round((index_max_cov-maxlag-1))/(Fs*1000)))
%     peak_value = max_cov;
%     peak_index = round((index_max_cov-maxlag-1))/(Fs*1000);
% elseif abs(max_cov) > abs(min_cov) & round((index_max_cov-maxlag-1))/(Fs*1000) < 0
%     text(-5000,.4,num2str(round((index_max_cov-maxlag-1))/(Fs*1000)))
%     peak_value = max_cov;
%     peak_index = round((index_max_cov-maxlag-1))/(Fs*1000);
% elseif abs(min_cov) > abs(max_cov) & round((index_min_cov-maxlag-1))/(Fs*1000) >= 0
%     text(5000,-.4,num2str(round((index_min_cov-maxlag-1))/(Fs*1000)))
%     peak_value = min_cov;
%     peak_index = round((index_min_cov-maxlag-1))/(Fs*1000);
% elseif abs(min_cov) > abs(max_cov) & round((index_min_cov-maxlag-1))/(Fs*1000) < 0
%     text(-5000,-.4,num2str(round((index_min_cov-maxlag-1))/(Fs*1000)))
%     peak_value = min_cov;
%     peak_index = round((index_min_cov-maxlag-1))/(Fs*1000);
% end

% assignin('base','peak_value',peak_value)
% assignin('base','peak_index',peak_index)

crosscorr.cov = cov;
crosscorr.lags = lags/(Fs*1000);
crosscorr.CI_bar_x = get(CI_bar,'XData');
crosscorr.CI_bar_y = get(CI_bar,'YData');


%%%%% Find Peaks
[a,b,c,d] = localpeak(crosscorr.cov',threshold,2);
% assignin('base','a',a)
% assignin('base','b',b)
% assignin('base','c',c)
% assignin('base','d',d)

max_covs = crosscorr.cov(b);
size_max_covs = size(max_covs,2);
max_lags = crosscorr.lags(b);
max_peaks(1:size_max_covs,1) = max_lags';
max_peaks(1:size_max_covs,2) = max_covs';
% assignin('base','max_covs',max_covs)
% assignin('base','max_lags',max_lags)

min_covs = crosscorr.cov(d);
size_min_covs = size(min_covs,2);
min_lags = crosscorr.lags(d);
min_peaks(1:size_min_covs,1) = min_lags';
min_peaks(1:size_min_covs,2) = min_covs';
% assignin('base','min_covs',min_covs)
% assignin('base','min_lags',min_lags)

% plot(max_lags,max_covs,'ro')
% plot(min_lags,min_covs,'go')





all_lags = [max_lags min_lags];
all_covs = [max_covs min_covs];

[sort_all_lags,all_lags_i] = sort(abs(all_lags));

assignin('base','all_lags',all_lags)
assignin('base','sort_all_lags',sort_all_lags)
assignin('base','all_lags_i',all_lags_i)

peaks_lags = all_lags(all_lags_i);
peaks_covs = all_covs(all_lags_i);
sorted_peaks(1:size(peaks_covs,2),1) = peaks_lags';
sorted_peaks(1:size(peaks_covs,2),2) = peaks_covs';


assignin('base','peaks_lags',peaks_lags)
assignin('base','peaks_covs',peaks_covs)

plot(peaks_lags,peaks_covs,'ro')

crosscorr.peaks.thresh = threshold;
crosscorr.peaks.max = max_peaks;
crosscorr.peaks.min = min_peaks;
crosscorr.peaks.sorted = sorted_peaks;






% crosscorr.peak_value = peak_value;
% crosscorr.peak_latency = peak_index;

%plot(crosscorr.peak_latency,crosscorr.peak_value,'ro')
%assignin('base','crosscorr',crosscorr)
