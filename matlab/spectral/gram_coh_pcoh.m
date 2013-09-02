function mt_hal_coherence = gram_coh_pcoh(filename,rang,TS1,TS1_descrip,TS2,TS2_descrip,TS3,TS3_descrip);

%,NW,K,pad,window,winstep,FS,rang);

% Multitaper Time-Frequency Coherence
% function [sp1,sp2,coh]=mtcoherence(TS1,TS2,NW,K,pad,window,winstep);
% TS1 : input time series 1
% TS2 : input time series 2
% NW = time bandwidth parameter (e.g. 3 or 4)
% K = number of data tapers kept, usually 2*NW -1 (e.g. 5 or 7 for above)
% pad = padding for individual window. Usually, choose power
% of two greater than but closest to size of moving window.
% window = length of moving window
% winstep = number of of timeframes between successive windows


% EEG_L = load_PCDX_SE('/F/Data/Raw/viv05/viv0510d.all','1-10',1);
% EEG_L_r = resample(EEG_L,400,10000);
% TS3 = EEG_L_r;
% TS3_descrip = 'EEG_L';
% 
% GC = load_PCDX_SE('/F/Data/Raw/viv05/viv0510d.all','1-10',3);
% GC_r = resample(GC,400,10000);
% TS2 = GC_r;
% TS2_descrip = 'GC';
% 
% DCN = load_PCDX_SE('/F/Data/Raw/viv05/viv0510d.all','1-10',4);
% spike_times = findspikes_win_SE(DCN,10,{-300 -75 .1 1},0);
% DCN_r = resamp_spike_times(spike_times,.4,1e4);
% TS1 = DCN_r;
% TS1_descrip = 'DCN';
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Multitapered Time-Frequency Coherence %%%%%%%%%%%%%%

%%% Set parameters
NW = 4;
K = 7;
pad = 4096;
window = 4096;
winstep = 512;
FS = 400;
%rang = [0 200];


TS1_mt = TS1(:)';
TS2_mt = TS2(:)';
TS3_mt = TS3(:)';
[E V] = dpss(window,NW,'calc');

[dum N]=size(TS1_mt);
sp1 = zeros(round((N-window)/winstep), pad/2);
sp2 = sp1;
sp3 = sp1;
coh12 = sp1;
coh13 = sp1;
coh23 = sp1;
for j = 1:((N-window)/winstep)
 
eJ = zeros(1,pad);

TSM1=TS1_mt((j-1)*winstep+[1:window])';
TSM1=TSM1-mean(TSM1);

TSM2=TS2_mt((j-1)*winstep+[1:window])';
TSM2=TSM2-mean(TSM2);

TSM3=TS3_mt((j-1)*winstep+[1:window])';
TSM3=TSM3-mean(TSM3);

J1=fft(TSM1(:,ones(1,K)).*(E(:,1:K)),pad)';
J2=fft(TSM2(:,ones(1,K)).*(E(:,1:K)),pad)';
J3=fft(TSM3(:,ones(1,K)).*(E(:,1:K)),pad)';


eJ1=real(sum(J1.*conj(J1),1));

eJ2=real(sum(J2.*conj(J2),1));

eJ3=real(sum(J3.*conj(J3),1));

eJ12=sum(J1.*conj(J2),1);
eJ13=sum(J1.*conj(J3),1);
eJ23=sum(J2.*conj(J3),1);


sp1(j,:)=eJ1(1:pad/2)/K;
sp2(j,:)=eJ2(1:pad/2)/K;
sp3(j,:)=eJ3(1:pad/2)/K;

% Calculate coherence
coh12(j,:)=eJ12(1:pad/2)./sqrt(eJ1(1:pad/2).*eJ2(1:pad/2));
coh13(j,:)=eJ13(1:pad/2)./sqrt(eJ1(1:pad/2).*eJ3(1:pad/2));
coh23(j,:)=eJ23(1:pad/2)./sqrt(eJ2(1:pad/2).*eJ3(1:pad/2));
% assignin('base','coh12',coh12)
% assignin('base','coh13',coh13)
% assignin('base','coh23',coh23)


end


% Calculate partial coherence
pcoh_12_p_3 = (coh12-coh13.*conj(coh23))./sqrt((1-abs(coh13).^2).*(1-abs(coh23).^2));
pcoh_12_p_3 = pcoh_12_p_3';
assignin('base','pcoh_12_p_3',pcoh_12_p_3)


pcoh_13_p_2 = (coh13-coh12.*coh23)./sqrt((1-abs(coh12).^2).*(1-abs(coh23).^2));
pcoh_13_p_2 = pcoh_13_p_2';
assignin('base','pcoh_13_p_2',pcoh_13_p_2)


pcoh_23_p_1 = (coh23-coh13.*conj(coh12))./sqrt((1-abs(coh12).^2).*(1-abs(coh13).^2));
pcoh_23_p_1 = pcoh_23_p_1';
assignin('base','pcoh_23_p_1',pcoh_23_p_1)


sp1 = sp1';
sp2 = sp2';
sp3 = sp3';
coh12 = coh12';
coh13 = coh13';
coh23 = coh23';

% if nargin > 8
len=size(sp1,2);
freq=FS*(0:pad/2-1)/pad;
assignin('base','freq',freq)
timebase=winstep/FS*(.5+(0:len-1));
assignin('base','timebase',timebase)
rang;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Halliday coherence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TS1 vs TS2 %%%%%
NFFT = 1024;
Fs = 400;
WINDOW = ones(NFFT,1)./sqrt(NFFT);

%%% Calculate coherence
for i = 1:size(TS1,2)
    temp= psd(TS1(:,i),NFFT,Fs,WINDOW,0,'mean');
    spctx(:,i)=temp;
    temp= psd(TS2(:,i),NFFT,Fs,WINDOW,0,'mean');
    spcty(:,i)=temp;
    [temp,F12] = csd(TS1(:,i),TS2(:,i),NFFT,Fs,WINDOW,0,'mean');
    cs(:,i)=temp;
end

mean_spctx = mean(spctx,2);
assignin('base','mean_spctx',mean_spctx)
mean_spcty = mean(spcty,2);
assignin('base','mean_spcty',mean_spcty)
mean_cs = mean(cs,2);

hal_R12 = mean_cs./sqrt(mean_spctx.*mean_spcty);
hal_R_phase12 = angle(hal_R12);
hal_coh12 = abs(hal_R12).^2;

% Calculate estimate of the upper 95% confidence limit
L = prod(size(TS1))/NFFT;
assignin('base','L',L)
CL12 = 1-(0.05).^(1/(L-1));


% Calculate partial coherence
Rxy12 = coherency_p(TS1,TS2,NFFT,Fs);
Rxz12 = coherency_p(TS1,TS3,NFFT,Fs);
Ryz12 = coherency_p(TS2,TS3,NFFT,Fs);
Rxy_z_12 = (Rxy12.R-Rxz12.R.*conj(Ryz12.R))./sqrt((1-abs(Rxz12.R).^2).*(1-abs(Ryz12.R).^2));
freq_pcoh12 = Rxy12.freq;

% Calculate estimate of the upper 95% confidence limit
L = prod(size(TS1))/NFFT;
assignin('base','L',L)
r = 1; 
CL_p_12 = 1-(0.05).^(1/(L-r-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TS1 vs TS3 %%%%%
NFFT = 1024;
Fs = 400;
WINDOW = ones(NFFT,1)./sqrt(NFFT);

for i = 1:size(TS1,2)
    temp= psd(TS1(:,i),NFFT,Fs,WINDOW,0,'mean');
    spctx(:,i)=temp;
    temp= psd(TS3(:,i),NFFT,Fs,WINDOW,0,'mean');
    spcty(:,i)=temp;
    [temp,F13] = csd(TS1(:,i),TS3(:,i),NFFT,Fs,WINDOW,0,'mean');
    cs(:,i)=temp;
end

mean_spctx = mean(spctx,2);
assignin('base','mean_spctx',mean_spctx)
mean_spcty = mean(spcty,2);
assignin('base','mean_spcty',mean_spcty)
mean_cs = mean(cs,2);

hal_R13 = mean_cs./sqrt(mean_spctx.*mean_spcty);
hal_R_phase13 = angle(hal_R13);
hal_coh13 = abs(hal_R13).^2;

% Calculate estimate of the upper 95% confidence limit
L = prod(size(TS1))/NFFT;
assignin('base','L',L)
CL13 = 1-(0.05).^(1/(L-1));


% Calculate partial coherence
Rxy13 = coherency_p(TS1,TS2,NFFT,Fs);
Rxz13 = coherency_p(TS1,TS3,NFFT,Fs);
Ryz13 = coherency_p(TS2,TS3,NFFT,Fs);
Rxz_y_13 = (Rxz13.R-Rxy13.R.*Ryz13.R)./sqrt((1-abs(Rxy13.R).^2).*(1-abs(Ryz13.R).^2));
freq_pcoh13 = Rxz13.freq;

% Calculate estimate of the upper 95% confidence limit
L = prod(size(TS1))/NFFT;
assignin('base','L',L)
r = 1; 
CL_p_13 = 1-(0.05).^(1/(L-r-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TS2 vs TS3 %%%%%
NFFT = 1024;
Fs = 400;
WINDOW = ones(NFFT,1)./sqrt(NFFT);

for i = 1:size(TS2,2)
    temp= psd(TS2(:,i),NFFT,Fs,WINDOW,0,'mean');
    spctx(:,i)=temp;
    temp= psd(TS3(:,i),NFFT,Fs,WINDOW,0,'mean');
    spcty(:,i)=temp;
    [temp,F23] = csd(TS2(:,i),TS3(:,i),NFFT,Fs,WINDOW,0,'mean');
    cs(:,i)=temp;
end

mean_spctx = mean(spctx,2);
assignin('base','mean_spctx',mean_spctx)
mean_spcty = mean(spcty,2);
assignin('base','mean_spcty',mean_spcty)
mean_cs = mean(cs,2);

hal_R23 = mean_cs./sqrt(mean_spctx.*mean_spcty);
hal_R_phase23 = angle(hal_R23);
hal_coh23 = abs(hal_R23).^2;

% Calculate estimate of the upper 95% confidence limit
L = prod(size(TS2))/NFFT;
assignin('base','L',L)
CL23 = 1-(0.05).^(1/(L-1));

pcoh_23_p_1 = (coh23-coh13.*conj(coh12))./sqrt((1-abs(coh12).^2).*(1-abs(coh13).^2));
% Calculate partial coherence
Rxy23 = coherency_p(TS1,TS2,NFFT,Fs);
Rxz23 = coherency_p(TS1,TS3,NFFT,Fs);
Ryz23 = coherency_p(TS2,TS3,NFFT,Fs);
Ryz_x_23 = (Ryz23.R-Rxz23.R.*conj(Rxy23.R))./sqrt((1-abs(Rxy23.R).^2).*(1-abs(Rxz23.R).^2));
freq_pcoh23 = Ryz23.freq;

% Calculate estimate of the upper 95% confidence limit
L = prod(size(TS1))/NFFT;
assignin('base','L',L)
r = 1; 
CL_p_23 = 1-(0.05).^(1/(L-r-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coh_01.hal_coh = coh;
% coh_01.hal_freq = F;
% plot(coh.freq,coh.coh); hold on
% plot(max_coh_range_freq,max_coh_range,'ro')
% x_lim = get(gca,'xlim');
% y_lim = get(gca,'ylim');
% plot([range(1) range(2)],[coh.CI coh.CI],'r');
% axis([range(1) range(2) 0 1]);
% ylabel('Coherence');
% xlabel('Frequency (Hz)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TS1_TS2_str_01 = ['mt_hal_coherence.',TS1_descrip,'_',TS2_descrip,'.mt_coh = abs(coh12).^2;'];
eval(TS1_TS2_str_01)
TS1_TS2_str_02 = ['mt_hal_coherence.',TS1_descrip,'_',TS2_descrip,'.mt_pcoh = abs(pcoh_12_p_3).^2;'];
eval(TS1_TS2_str_02)
TS1_TS2_str_03 = ['mt_hal_coherence.',TS1_descrip,'_',TS2_descrip,'.mt_timebase = timebase;'];
eval(TS1_TS2_str_03)
TS1_TS2_str_04 = ['mt_hal_coherence.',TS1_descrip,'_',TS2_descrip,'.mt_freq = freq;'];
eval(TS1_TS2_str_04)
TS1_TS2_str_05 = ['mt_hal_coherence.',TS1_descrip,'_',TS2_descrip,'.hal_ccoh = hal_coh12;'];
eval(TS1_TS2_str_05)
TS1_TS2_str_06 = ['mt_hal_coherence.',TS1_descrip,'_',TS2_descrip,'.hal_Rphase = hal_R_phase12;'];
eval(TS1_TS2_str_06)
TS1_TS2_str_07 = ['mt_hal_coherence.',TS1_descrip,'_',TS2_descrip,'.hal_cfreq = F12;'];
eval(TS1_TS2_str_07)
TS1_TS2_str_08 = ['mt_hal_coherence.',TS1_descrip,'_',TS2_descrip,'.hal_cCI = CL12;'];
eval(TS1_TS2_str_08)
TS1_TS2_str_09 = ['mt_hal_coherence.',TS1_descrip,'_',TS2_descrip,'.hal_pcoh = abs(Rxy_z_12).^2;'];
eval(TS1_TS2_str_09)
TS1_TS2_str_10 = ['mt_hal_coherence.',TS1_descrip,'_',TS2_descrip,'.hal_pfreq = freq_pcoh12;'];
eval(TS1_TS2_str_10)
TS1_TS2_str_11 = ['mt_hal_coherence.',TS1_descrip,'_',TS2_descrip,'.hal_pCI = CL_p_12;'];
eval(TS1_TS2_str_11)

TS1_TS3_str_01 = ['mt_hal_coherence.',TS1_descrip,'_',TS3_descrip,'.mt_coh = abs(coh13).^2;'];
eval(TS1_TS3_str_01)
TS1_TS3_str_02 = ['mt_hal_coherence.',TS1_descrip,'_',TS3_descrip,'.mt_pcoh = abs(pcoh_13_p_2).^2;'];
eval(TS1_TS3_str_02)
TS1_TS3_str_03 = ['mt_hal_coherence.',TS1_descrip,'_',TS3_descrip,'.mt_timebase = timebase;'];
eval(TS1_TS3_str_03)
TS1_TS3_str_04 = ['mt_hal_coherence.',TS1_descrip,'_',TS3_descrip,'.mt_freq = freq;'];
eval(TS1_TS3_str_04)
TS1_TS3_str_05 = ['mt_hal_coherence.',TS1_descrip,'_',TS3_descrip,'.hal_ccoh = hal_coh13;'];
eval(TS1_TS3_str_05)
TS1_TS3_str_06 = ['mt_hal_coherence.',TS1_descrip,'_',TS3_descrip,'.hal_Rphase = hal_R_phase13;'];
eval(TS1_TS3_str_06)
TS1_TS3_str_07 = ['mt_hal_coherence.',TS1_descrip,'_',TS3_descrip,'.hal_cfreq = F13;'];
eval(TS1_TS3_str_07)
TS1_TS3_str_08 = ['mt_hal_coherence.',TS1_descrip,'_',TS3_descrip,'.hal_cCI = CL13;'];
eval(TS1_TS3_str_08)
TS1_TS3_str_09 = ['mt_hal_coherence.',TS1_descrip,'_',TS3_descrip,'.hal_pcoh = abs(Rxz_y_13).^2;'];
eval(TS1_TS3_str_09)
TS1_TS3_str_10 = ['mt_hal_coherence.',TS1_descrip,'_',TS3_descrip,'.hal_pfreq = freq_pcoh13;'];
eval(TS1_TS3_str_10)
TS1_TS3_str_11 = ['mt_hal_coherence.',TS1_descrip,'_',TS3_descrip,'.hal_pCI = CL_p_13;'];
eval(TS1_TS3_str_11)

TS2_TS3_str_01 = ['mt_hal_coherence.',TS2_descrip,'_',TS3_descrip,'.mt_coh = abs(coh23).^2;'];
eval(TS2_TS3_str_01)
TS2_TS3_str_02 = ['mt_hal_coherence.',TS2_descrip,'_',TS3_descrip,'.mt_pcoh = abs(pcoh_23_p_1).^2;'];
eval(TS2_TS3_str_02)
TS2_TS3_str_03 = ['mt_hal_coherence.',TS2_descrip,'_',TS3_descrip,'.mt_timebase = timebase;'];
eval(TS2_TS3_str_03)
TS2_TS3_str_04 = ['mt_hal_coherence.',TS2_descrip,'_',TS3_descrip,'.mt_freq = freq;'];
eval(TS2_TS3_str_04)
TS2_TS3_str_05 = ['mt_hal_coherence.',TS2_descrip,'_',TS3_descrip,'.hal_ccoh = hal_coh23;'];
eval(TS2_TS3_str_05)
TS2_TS3_str_06 = ['mt_hal_coherence.',TS2_descrip,'_',TS3_descrip,'.hal_Rphase = hal_R_phase23;'];
eval(TS2_TS3_str_06)
TS2_TS3_str_07 = ['mt_hal_coherence.',TS2_descrip,'_',TS3_descrip,'.hal_cfreq = F23;'];
eval(TS2_TS3_str_07)
TS2_TS3_str_08 = ['mt_hal_coherence.',TS2_descrip,'_',TS3_descrip,'.hal_cCI = CL23;'];
eval(TS2_TS3_str_08)
TS2_TS3_str_09 = ['mt_hal_coherence.',TS2_descrip,'_',TS3_descrip,'.hal_pcoh = abs(Ryz_x_23).^2;'];
eval(TS2_TS3_str_09)
TS2_TS3_str_10 = ['mt_hal_coherence.',TS2_descrip,'_',TS3_descrip,'.hal_pfreq = freq_pcoh23;'];
eval(TS2_TS3_str_10)
TS2_TS3_str_11 = ['mt_hal_coherence.',TS2_descrip,'_',TS3_descrip,'.hal_pCI = CL_p_23;'];
eval(TS2_TS3_str_11)

%assignin('base','mt_hal_coherence',mt_hal_coherence)





figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Plot Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%% Signal 2 vs Signal 3
   subplot(4,15,1:5)
   imagesc(timebase,freq,abs(coh23).^2)
   axis xy
   set(gca,'ylim',rang,'xticklabel',[])
   title([filename,' ',TS2_descrip,' vs ',TS3_descrip],'FontSize',9)
   ylabel('freq (Hz)')
   subplot(4,15,6:7)
   plot(hal_coh23,F23); hold on
   plot([CL23 CL23],[rang(1) rang(2)],'r')
   %plot(abs(mean(coh23,2)).^2,freq,'b')
   %set(gca,'ylim',rang,'xlim',[0 1],'xticklabel',[])
   set(gca,'ylim',rang,'xlim',[0 1])
   
   subplot(4,15,8:13)
   imagesc(timebase,freq,abs(pcoh_23_p_1).^2)
   axis xy
   title([TS2_descrip,' vs ',TS3_descrip,' with ',TS1_descrip,' as predictor'],'FontSize',9)
   set(gca,'yticklabel',[],'ylim',rang,'xticklabel',[],'Position',[0.55 0.767258 0.248 0.157742])
   subplot(4,15,14:15)
   plot(abs(Ryz_x_23).^2,freq_pcoh23); hold on
   plot([CL_p_23 CL_p_23],[rang(1) rang(2)],'r')
   %plot(abs(mean(pcoh_23_p_1,2)).^2,freq,'b')
   %set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1],'xticklabel',[])
   set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1])
   
   %%%%% Signal 1 vs Signal 3
   subplot(4,15,16:20)
   imagesc(timebase,freq,abs(coh13).^2)
   axis xy
   set(gca,'ylim',rang,'xticklabel',[],'Position',[0.13 0.548172 0.248461 0.157742])
   title([TS1_descrip,' vs ',TS3_descrip],'FontSize',9)
   ylabel('freq (Hz)')
   subplot(4,15,21:22)
   plot(hal_coh13,F13); hold on
   plot([CL13 CL13],[rang(1) rang(2)],'r')
   %plot(abs(mean(coh13,2)).^2,freq,'b')
   %set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1],'xticklabel',[])
   set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1])
   
   subplot(4,15,23:28)
   imagesc(timebase,freq,abs(pcoh_13_p_2).^2)
   axis xy
   title([TS1_descrip,' vs ',TS3_descrip,' with ',TS2_descrip,' as predictor'],'FontSize',9)
   set(gca,'yticklabel',[],'xticklabel',[],'ylim',rang,'Position',[0.55 0.548172 0.248 0.157742])
   subplot(4,15,29:30)
   plot(abs(Rxz_y_13).^2,freq_pcoh13); hold on
   plot([CL_p_13 CL_p_13],[rang(1) rang(2)],'r')
   %plot(abs(mean(pcoh_13_p_2,2)).^2,freq,'b')
   %set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1],'xticklabel',[])
   set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1])
   
   %%%%% Signal 1 vs Signal 2
   subplot(4,15,31:35)
   imagesc(timebase,freq,abs(coh12).^2)
   axis xy
   set(gca,'ylim',rang,'Position',[0.13 0.329 0.248461 0.157742])
   title([TS1_descrip,' vs ',TS2_descrip],'FontSize',9)
   ylabel('freq (Hz)')
   xlabel('time (s)')
   subplot(4,15,36:37)
   plot(hal_coh12,F12); hold on
   plot([CL12 CL12],[rang(1) rang(2)],'r')
   %plot(abs(mean(coh12,2)).^2,freq,'b')
   xlabel('Coherence','FontSize',7)
   set(gca,'ylim',rang,'xlim',[0 1],'Position',[0.392534 0.329 0.0855836 0.157742])
   
   subplot(4,15,38:43)
   imagesc(timebase,freq,abs(pcoh_12_p_3).^2)
   axis xy
   title([TS1_descrip,' vs ',TS2_descrip,' with ',TS3_descrip,' as predictor'],'FontSize',9)
   set(gca,'yticklabel',[],'ylim',rang,'Position',[0.55 0.329 0.248 0.157742])
   xlabel('time (s)')
   subplot(4,15,44:45)
   plot(abs(Rxy_z_12).^2,freq_pcoh12); hold on
   plot([CL_p_12 CL_p_12],[rang(1) rang(2)],'r')
   %plot(abs(mean(pcoh_12_p_3,2)).^2,freq,'b')
   xlabel(['Partial  ';'Coherence'],'FontSize',7)
   set(gca,'ylim',rang,'xlim',[0 1],'Position',[0.812588 0.329 0.0924119 0.157742])
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

