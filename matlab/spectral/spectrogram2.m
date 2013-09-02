function spec_gram = spectrogram_SE(signal,rang,winwidth,threshold,pos,sp1a,sp2a,sp3a)

% spectrogram This function calculates and plots a time-resolved
% power spectrum of the given signal.  The function will return 
% frequencies up to the Nyquist frequency, but will only plot
% the range of frequencies entered (e.g., [0 30]).  A mean
% power spectrum will also be plotted over all times. Enter
% a window width for smoothing the power spectrum.  Also, enter
% a threshold for finding peaks in the mean power spectrum.  To plot:
%     - 
%     - 
% 
% spec_gram = spectrogram(signal,range,winwidth,threshold,pos,sp1a,sp2a,sp3a)

TS1 = signal;
NW = 4;
K = 7;
pad = 4096;
window = 4096;
winstep = 512;
FS = 400;
%rang = [0,200];



TS1 = TS1(:)';
[E V] = dpss(window,NW,'calc');
%assignin('base','E',E)
%assignin('base','V',V)
[dum N] = size(TS1);
sp1 = zeros(round((N-window)/winstep), pad/2);

for j = 1:((N-window)/winstep)
    eJ = zeros(1,pad);
    TSM1 = TS1((j-1)*winstep+[1:window])';
    TSM1 = TSM1-mean(TSM1);
    J1 = fft(TSM1(:,ones(1,K)).*(E(:,1:K)),pad)';
    eJ1 = real(sum(J1.*conj(J1),1));
    sp1(j,:) = eJ1(1:pad/2)/K;
end

sp1 = sp1';
assignin('base','sp1',sp1)

len = size(sp1,2);
freq = FS*(0:pad/2-1)/pad;
assignin('base','freq',freq)
assignin('base','sp1',sp1)
timebase = winstep/FS*(.5+(0:len-1)); 
assignin('base','timebase',timebase)
rang;

% [max_ps_y,max_ps_i] = max(mean(10*log10(sp1),2));
% assignin('base','max_ps_y',max_ps_y)
% assignin('base','max_ps_i',max_ps_i)
% max_freq = freq(max_ps_i);
% assignin('base','max_freq',max_freq)

% Smooth power spectrum
power_spectrum = mean(10*log10(sp1),2);
power_spectrum_s = smooth(power_spectrum,winwidth);

% Find peaks
%power_spectrum = mean(10*log10(sp1),2);
find_freq_rang = find(freq >= 0 & freq <= 30);
freq_rang = freq(find_freq_rang);
power_spectrum_rang = power_spectrum_s(find_freq_rang);
[a,b,c,d] = localpeak(power_spectrum_rang,threshold,2);
max_freq_rang = freq_rang(b)';
max_power_spectrum_rang = a;
peaks_rang(1:size(max_freq_rang,1),1) = max_freq_rang;
peaks_rang(1:size(max_freq_rang,1),2) = max_power_spectrum_rang;
% assignin('base','a',a)
% assignin('base','b',b)
% assignin('base','c',c)
% assignin('base','d',d)
% assignin('base','freq_rang',freq_rang)
% assignin('base','power_spectrum_rang',power_spectrum_rang)
% assignin('base','max_freq_rang',max_freq_rang)
% assignin('base','max_power_spectrum_rang',max_power_spectrum_rang)
% assignin('base','peaks',peaks)



% Assign outputs
spec_gram.spec = sp1;
spec_gram.freq = freq;
spec_gram.timebase = timebase;
spec_gram.peaks = peaks_rang;

% plot
if nargin == 3
    figure_posi = [(pos-1)*5+1:pos*5];
    subplot(13,15,figure_posi(1:3))
    imagesc(timebase,freq,10*log10(sp1))
    axis xy
    set(gca,'ylim',rang)
    ylabel('freq (Hz)')
    xlabel('time (s)')
    subplot(13,15,figure_posi(4:5))
    plot(mean(10*log10(sp1),2),freq); hold on
    plot(max_ps_y,max_freq,'ro')
    set(gca,'FontSize',7,'XTick',[],'YTick',[])
    set(gca,'yticklabel',[],'ylim',rang,'FontSize',7)
    subplot(13,15,figure_posi(1:3))
elseif nargin == 8
    divide_columns = sp2a*5;
    recalc_fig_no = sp3a*5;
    figure_pos = [recalc_fig_no-4:recalc_fig_no];
    subplot(sp1a,divide_columns,figure_pos(1:3))
    imagesc(timebase,freq,10*log10(sp1))
    axis xy
    set(gca,'ylim',rang)
    ylabel('freq (Hz)')
    xlabel('time (s)')
    subplot(sp1a,divide_columns,figure_pos(4:5))
    %plot(mean(10*log10(sp1),2),freq); hold on
    plot(power_spectrum_s,freq); hold on
    plot(spec_gram.peaks(:,2),spec_gram.peaks(:,1),'ro')
    %plot(max_ps_y,max_freq,'ro')
    %set(gca,'FontSize',7,'XTick',[],'YTick',[])
    set(gca,'yticklabel',[],'ylim',rang,'FontSize',7)
    subplot(sp1a,divide_columns,figure_pos(1:3))
else
    subplot(4,15,1:5)
    imagesc(timebase,freq,10*log10(sp1))
    axis xy
    set(gca,'ylim',rang)
    ylabel('freq (Hz)')
    xlabel('time (s)')
    subplot(4,15,6:7)
    plot(mean(10*log10(sp1),2),freq); hold on
    plot(max_ps_y,max_freq,'ro')
    set(gca,'yticklabel',[],'ylim',rang)
end