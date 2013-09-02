function p_coh = partial_coherency_new(X,Y,Z,NFFT,Fs,range,plotit,sp1,sp2,sp3);

% partial_coherency_new This function calculates the partial coherency between
% signals X and Y with signal Z serving as the predictor.  Partial coherency 
% tests the hypothesis that the coherence between X and Y is entirely due to 
% the common influence of Z, in which case the partial coherence will be zero 
% at a particular frequency band.  If the coherence between X and Y is entirely 
% independent of Z, then the partial coherence will be 1 at a particular frequency 
% (Halliday et al, 1995). The signals should be first resampled to 400Hz. Specify 
% the frequency range you wish the algorithm to analyze as [min_freq max_freq].  
% NOVERLAP has been set at 0 and FLAG has been set to 'mean'. Enter '1' to plot 
% the resulting graph; otherwise enter '0'.  If plot is set to '1', you may also 
% specify subplot coordinates
%
% p_coh = partial_coherency(X,Y,Z,NFFT,Fs,plotit,sp1,sp2,sp3);
%
% Example 1: EEG_L = load_PCDX_SE('/Raw/viv05/viv0511a.all','1',1);
%            EEG_L_r = resample(EEG_L,400,10000);
%            EEG_R = load_PCDX_SE('/Raw/viv05/viv0511a.all','1',2);
%            EEG_R_r = resample(EEG_R,400,10000);
%            GC = load_PCDX_SE('/Raw/viv05/viv0511a.all','1',3);
%            GC_r = resample(GC,400,10000);
%            p_coh = partial_coherency_new(EEG_L_r,EEG_R_r,GC_r,1024,400,[0 25],1,2,2,1);



% % Resample signals at 400Hz
% X = resample(X,400,Fs);
% %Y = resample(Y,400,Fs);
% Z = resample(Z,400,Fs);

% Calculate partial coherence
Rxy = coherency_p(X,Y,NFFT,Fs);
Rxz = coherency_p(X,Z,NFFT,Fs);
Ryz = coherency_p(Y,Z,NFFT,Fs);
Rxy_z = (Rxy.R-Rxz.R.*conj(Ryz.R))./sqrt((1-abs(Rxz.R).^2).*(1-abs(Ryz.R).^2));

% Calculate estimate of the upper 95% confidence limit
L = prod(size(X))/NFFT;
assignin('base','L',L)
r = 1; % r = number of predictors, which is 1
CL = 1-(0.05).^(1/(L-r-1));


% % Calculate estimate of the upper 95% confidence limit
% L = prod(size(X))/NFFT;
% assignin('base','L',L)
% CL = 1-(0.05).^(1/(L-1));



% Send appropriate values to structure
p_coh.freq = Rxy.freq;
p_coh.coh = Rxy.coh;
p_coh.p_coh = abs(Rxy_z).^2;
p_coh.cl = CL;
% p_coh.Rxz = Rxz.coh;
% p_coh.Ryz = Ryz.coh;

%assignin('base','F',F)
find_F_range = find(p_coh.freq <= range(2));
%assignin('base','find_F_range',find_F_range)
coh_range = p_coh.coh(find_F_range);
%assignin('base','coh_range',coh_range)
[max_coh_range,max_coh_i_range] = max(coh_range);
%assignin('base','max_coh_i_range',max_coh_i_range)
%assignin('base','max_coh_range',max_coh_range)
max_coh_range_freq = p_coh.freq(max_coh_i_range);
%assignin('base','max_coh_range_freq',max_coh_range_freq)
max_p_coh_range = p_coh.p_coh(max_coh_i_range);
assignin('base','max_p_coh_range',max_p_coh_range)
p_coh.peak_coh_range_x = max_coh_range_freq;
p_coh.peak_coh_range_y = max_coh_range;
p_coh.peak_p_coh_range_y = max_p_coh_range;
coh_range_delta = max_coh_range - max_p_coh_range;
p_coh.peak_coh_range_delta = coh_range_delta;

% Plot
  
if plotit == 1
    if nargin == 10
        subplot(sp1,sp2,sp3)
    end
%     figure
%     subplot(4,1,3)
    
    coh_plot = plot(p_coh.freq,p_coh.coh,'b','LineWidth',2); hold on
    p_coh_plot = plot(p_coh.freq,p_coh.p_coh,'r')
%     leg = legend('Coherency X|Y','Coherency X|Y w/ Z predictor')
%     set(leg,'FontSize',7)
    plot(max_coh_range_freq,max_coh_range,'bo','Linewidth',2)
    plot(max_coh_range_freq,max_p_coh_range,'ro','Linewidth',2)
    x_lim = get(gca,'xlim');
    y_lim = get(gca,'ylim');
    plot([x_lim(1) x_lim(2)],[p_coh.cl p_coh.cl],'r','LineWidth',2)
    axis([range(1) range(2) 0 1])
    ylabel('Partial Coherence');
    xlabel('Frequency (Hz)');
end









