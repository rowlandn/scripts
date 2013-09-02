function coh = coherency3(X,Y,NFFT,Fs,range,plotit,sp1,sp2,sp3);

% coherency This function calculates the coherence between two signals.
% Coherence is a measure of the linear association of two signals in 
% the frequency domain (Halliday et al, 1995).  Coherence is based on 
% the power spectra of the two signals and is bounded between 0 and 1, 
% with 0 signifying no relationship between the signals and 1 signifying 
% that the two signals are identical.  The two signals must be represented 
% as column vectors.  The signals should be resampled to 400Hz using the original 
% sampling rate (Fs) of the recordings in Hz. In this case, use NFFT = 1024.
% NOVERLAP has been set at 0 and FLAG has been set to 'mean'.  Enter a range
% of frequencies in Hz to plot, e.g., [0 25] (keep in mind, though, that the 
% function will return coherences up to the Nyquist frequency. Enter '1' to 
% plot the resulting graph; otherwise enter '0'.  If plotit is 1, the peak
% of the graph in the frequency range entered will be identified with a red
% circle.
%
% coh = coherency3(X,Y,NFFT,Fs,range,plotit,sp1,sp2,sp3);
% 
% The structure coh will return:
%                              coh: coherence (y) values
%                             freq: frequency (x) values
%                               CI: upper 95% confidence limit
%                                R: ?
%                          R_phase: ?
%                       mean_spctx: power spectrum of signal X
%                       mean_spcty: power spectrum of signal Y
%          max_coh_range_above_sig: peak coherence above 
%                                   significance in frequency 
%                                   range specified
%     max_coh_range_freq_above_sig: frequency corresponding to
%                                   peak coherence above significance 
%                                   in frequency range specified
% 
% Example 1: 
% EEG_L = load_PCDX_SE('R:\viv05\viv0511a.all','1',1);
%            EEG_L_r = resample(EEG_L,400,10000);
%            GC = load_PCDX_SE('R:\viv05\viv0511a.all','1',3);
%            GC_r = resample(GC,400,10000);
%            coh = coherency3(EEG_L_r,GC_r,1024,400,[0 25],1);
%
% Example 2: EEG_L = load_PCDX_SE('/F/Data/Raw/viv05/viv0511a.all','1',1);
%            EEG_L_r = resample(EEG_L,400,10000);
%            GC = load_PCDX_SE('/F/Data/Raw/viv05/viv0511a.all','1',3);
%            GC_r = resample(GC,400,10000);
%            coh = coherency3(EEG_L_r,GC_r,1024,400,[0 25],1,2,2,1);
%             (if you want to specify a certain subplot,
%                       use this notation)




% % Resample signals at 400Hz
% X = resample(X,400,OFs);
% Y = resample(Y,400,OFs);

% Calculate coherence
WINDOW = ones(NFFT,1)./sqrt(NFFT);

for i = 1:size(X,2)

temp= psd(X(:,i),NFFT,Fs,WINDOW,0,'mean');
spctx(:,i)=temp;
temp= psd(Y(:,i),NFFT,Fs,WINDOW,0,'mean');
spcty(:,i)=temp;
[temp,F] = csd(X(:,i),Y(:,i),NFFT,Fs,WINDOW,0,'mean');
cs(:,i)=temp;
end

mean_spctx = mean(spctx,2);
assignin('base','mean_spctx',mean_spctx)
mean_spcty = mean(spcty,2);
assignin('base','mean_spcty',mean_spcty)
mean_cs = mean(cs,2);

R = mean_cs./sqrt(mean_spctx.*mean_spcty);
coh = abs(R).^2;



% Calculate estimate of the upper 95% confidence limit
L = prod(size(X))/NFFT;
assignin('base','L',L)
CL = 1-(0.05).^(1/(L-1));

% Send appropriate values to structure

%assignin('base','F',F)
find_F_range = find(F <= range(2));
%assignin('base','find_F_range',find_F_range)
coh_range = coh(find_F_range);
%assignin('base','coh_range',coh_range)
[max_coh_range,max_coh_i_range] = max(coh_range);
%assignin('base','max_coh_i_range',max_coh_i_range)
%assignin('base','max_coh_range',max_coh_range)
max_coh_range_freq = F(max_coh_i_range);
%assignin('base','max_coh_range_freq',max_coh_range_freq)



coh.coh = coh;
coh.freq = F;
coh.CI = CL;
coh.R = R;
coh.R_phase = angle(R);
coh.mean_spctx = mean_spctx;
coh.mean_spcty = mean_spcty;

if max_coh_range >= CL & max_coh_range_freq <= range(2)
    coh.max_coh_range_above_sig = max_coh_range-CL;
    coh.max_coh_range_freq_above_sig = max_coh_range_freq;
else
    coh.max_coh_range_above_sig = 0;
    coh.max_coh_range_freq_above_sig = 0;
end




% Plot
if plotit == 1
    if nargin == 9
        subplot(sp1,sp2,sp3)
    else
    end
%     %figure
%     subplot(4,1,1)
%     semilogy(coh.freq,coh.mean_spctx)
%     ylabel([sig_X_descrip,' spectrum'])
%     title(filename)
%     
%     subplot(4,1,2)
%     semilogy(coh.freq,coh.mean_spcty)
%     ylabel([sig_Y_descrip,' spectrum'])
    
    %subplot(4,1,3)
    plot(coh.freq,coh.coh); hold on
    %plot(max_coh_range_freq,max_coh_range,'ro')
    x_lim = get(gca,'xlim');
    y_lim = get(gca,'ylim');
    %plot([range(1) range(2)],[coh.CI coh.CI],'r');
    axis([range(1) range(2) -3.5 3.5]);
    ylabel('Coherence');
    xlabel('Frequency (Hz)')
    
%     subplot(4,1,4)
%     plot(coh.freq,coh.R_phase)
%     ylabel('Coherency Phase')
%     xlabel('Frequency (Hz)');
end