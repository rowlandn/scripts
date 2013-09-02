% this script plots the frequency spectrum of files containing digitized
% ecog data.  The script looks for paired variables named "ecogrest" and
% "ecogactive" (recordings from same contact montage, during rest and active self-cued movement)
% and plots the time domain and frequency domain for each, side by side

% define variables
FREQ_LO = [8 33];
FREQ_HI = [76 101];
YLIM = [-200 200];
Fs=1000;        %FS = samples per unit time, in this case digitization at 1000 hz rate
%Gain=5097.7    %gain is now pulled directly pulled from patient file using
                %apmconv5, thus there is no need to hardwire gain
BIN = 1;        % bin used for FFT downsampling
SPEC_MAX = 100; % max end of spectral frequency for plotting


%% ecogrest time pre-downsample
z=length(ecogrest);
lastnum=(z-1)/Fs;
ecogresttime=0:1/Fs:lastnum;

subplot(4,2,1);
plot(ecogresttime,ecogrest); %plots the resting ecog in time domain, scale in microvolts
xlabel('time (seconds)');
ylabel('voltage (uV)');
ylim(YLIM);
title('ecog at rest in time domain')

Y=fft(ecogrest);
Y(1)=0; %for some reason the first number in the fft vector is the sum of all others, 
% and to avoid confusion in the plotting, we set that first number to zero
n=length(Y); 
P=(abs(Y).^2)/n; %calculates the power of the fft, which is square of its magnitude
P=P./max(P); %normalizes the power to the maximum power of 1 to facilate plotting
P1=P(1:floor(n./2));%pulls out only the first half of the spectrum up to nyquist frequency




n1=length(P1);
f=(0:n1-1)*(Fs/n); %sets the frequency or x scale to be correct range to correspond to fft
subplot(4,2,3);
plot(f,P1); %plots the power versus the frequency
axis([0 SPEC_MAX 0 1]); %expands the first 100 Hz of the spectrum
xlabel('frequency (Hz)');
ylabel('normalized power');
title('ecog frequencies at rest, 0-100 hz')

% next few lines calculates the percentage of total power in the low frequency (8-32 hz) and high
% frequency (76-100 hz) ends of the spectrum.  variables are named LFPr
% (low frequency power at rest) and HFPr (high frequency power at rest)

TP=sum(P1); %calculates total power in spectrum
Ilow=(f>FREQ_LO(1)) & (f<FREQ_LO(2)); % Ilow is a logical array pointing to low frequency range
LFPr=((sum(P1(Ilow)))/TP)*100; %sums all values of fft vector in low frequency range, normalizes to total power
Ihigh=(f>FREQ_HI(1)) & (f<FREQ_HI(2)); %Ihigh is a logical array pointing to high frequency range
HFPr=((sum(P1(Ihigh)))/TP)*100; %sums all values of fft vector in high frequency range, normalizes to total power
ratiorest=HFPr/LFPr;

%% ecog active pre-downsample
za=length(ecogactive);
lastnuma=(za-1)/Fs;
ecogactivetime=0:1/Fs:lastnuma;

subplot(4,2,2);
plot(ecogactivetime,ecogactive); %plots the active movement ecog in time domain
xlabel('time (seconds)');
ylabel('voltage (uV)');
ylim(YLIM);
title('ecog during movement in time domain')

Ya=fft(ecogactive);
Ya(1)=0; %for some reason the first number in the fft vector is the sum of all others, 
% and to avoid confusion in the plotting, we set that first number to zero
na=length(Ya); 
Pa=(abs(Ya).^2)/na; %calculates the power of the fft, which is square of its magnitude
Pa=Pa./max(Pa); %normalizes the power to the maximum power of 1 to facilate plotting
P1a=Pa(1:floor(na./2));%pulls out only the first half of the spectrum up to nyquist frequency
n1a=length(P1a);
fa=(0:n1a-1)*(Fs/na); %sets the frequency or x scale to be correct range to correspond to fft
subplot(4,2,4);
plot(fa,P1a); %plots the power versus the frequency
axis([0 SPEC_MAX 0 1]); %expands the first 100 Hz of the spectrum
xlabel('frequency (Hz)');
ylabel('normalized power');
title('ecog frequency during movement, 0-100 hz')


% next few lines calculates the percentage of total power in the low frequency (8-32 hz) and high
% frequency (76-100 hz) ends of the spectrum.  variables are named LFPr
% (low frequency power at rest) and HFPr (high frequency power at rest)

TPa=sum(P1a); %calculates total power in spectrum
Ilowa=(fa>FREQ_LO(1)) & (fa<FREQ_LO(2)); % Ilowa is a logical array pointing to low frequency range
LFPa=((sum(P1a(Ilowa)))/TPa)*100; %sums all values of fft vector in low frequency range, normalizes to total power
Ihigha=(fa>FREQ_HI(1)) & (fa<FREQ_HI(2)); %Ihigha is a logical array pointing to high frequency range
HFPa=((sum(P1a(Ihigha)))/TPa)*100; %sums all values of fft vector in high frequency range, normalizes to total power
ratioactive=HFPa/LFPa;

%% ecog rest post-downsample
% initialize
P2 = zeros(1,SPEC_MAX/BIN);

for i = 1:length(P2)
    I_bin = (f>=BIN*(i-1)) & (f<BIN*i);
    P2(i) = mean(P1(I_bin));
end
f2 = 0:BIN:SPEC_MAX-BIN;
subplot(4,2,5);
plot(f2,P2);
subplot(4,2,7);
semilogy(f2,P2);


%% ecog active post-downsample
% initialize
P2a = zeros(1,SPEC_MAX/BIN);

for i = 1:length(P2a)
    I_bin = (fa>=BIN*(i-1)) & (fa<BIN*i);
    P2a(i) = mean(P1a(I_bin));
end
f2a = 0:BIN:SPEC_MAX-BIN;
subplot(4,2,6);
plot(f2a,P2a);
subplot(4,2,8);
semilogy(f2a,P2a);
%% new plots
figure;
plot(f2,P2,'-b',...
    'LineWidth',2);
hold on
plot(f2a,P2a,'-r',...
    'LineWidth',2);
set(gca,'XLim',[0 100]);
ylm = get(gca,'YLim');
x_lo = [FREQ_LO(1) FREQ_LO(2) FREQ_LO(2) FREQ_LO(1)];
x_hi = [FREQ_HI(1) FREQ_HI(2) FREQ_HI(2) FREQ_HI(1)];
y = [ylm(1) ylm(1) ylm(2) ylm(2)];
fill(x_lo,y,'g',...
    'EdgeColor','none',...
    'FaceAlpha',0.3);
fill(x_hi,y,'y',...
    'EdgeColor','none',...
    'FaceAlpha',0.3);
hold off
legend('Rest','Hand Movement');

figure;
semilogy(f2,P2,'-b',...
    'LineWidth',2);
hold on
semilogy(f2a,P2a,'-r',...
    'LineWidth',2);
set(gca,'XLim',[0 100]);
ylm_log = get(gca, 'YLim');
y_log = [ylm_log(1) ylm_log(1) ylm_log(2) ylm_log(2)];
fill(x_lo,y_log,'g',...
    'EdgeColor','none',...
    'FaceAlpha',0.3);
fill(x_hi,y_log,'y',...
    'EdgeColor','none',...
    'FaceAlpha',0.3);
hold off
legend('Rest','Hand Movement');
