function coherogram4(filename,TS1,TS1_descrip,TS2,TS2_descrip,TS3,TS3_descrip);

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
NW = 4;
K = 7;
pad = 4096;
window = 4096;
winstep = 512;
FS = 400;
rang = [0,25];


TS1 = TS1(:)';
TS2 = TS2(:)';
TS3 = TS3(:)';
[E V] = dpss(window,NW,'calc');
%assignin('base','E',E)
%assignin('base','V',V)
[dum N]=size(TS1);
sp1 = zeros(round((N-window)/winstep), pad/2);
sp2 = sp1;
sp3 = sp1;
coh12 = sp1;
coh13 = sp1;
coh23 = sp1;
for j = 1:((N-window)/winstep)
 
eJ = zeros(1,pad);

TSM1=TS1((j-1)*winstep+[1:window])';
TSM1=TSM1-mean(TSM1);

TSM2=TS2((j-1)*winstep+[1:window])';
TSM2=TSM2-mean(TSM2);

TSM3=TS3((j-1)*winstep+[1:window])';
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
assignin('base','coh12',coh12)
assignin('base','coh13',coh13)
assignin('base','coh23',coh23)


end


% Calculate partial coherence
pcoh_12_p_3 = (coh12-coh13.*conj(coh23))./sqrt((1-abs(coh13).^2).*(1-abs(coh23).^2));
pcoh_12_p_3 = pcoh_12_p_3';
assignin('base','pcoh_12_p_3',pcoh_12_p_3)

% pcoh_13_p_2 = (coh13-coh12.*conj(coh23))./sqrt((1-abs(coh12).^2).*(1-abs(coh23).^2));
% pcoh_13_p_2 = pcoh_13_p_2';
% assignin('base','pcoh_13_p_2',pcoh_13_p_2)

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

%figure
%set(gcf,'Position',[20 50 1000 500])

   %%%%% Signal 1 vs Signal 2
   subplot(4,15,1:5)
   imagesc(timebase,freq,abs(coh12).^2)
   axis xy
   set(gca,'ylim',rang,'xticklabel',[])
   title([TS1_descrip,' vs ',TS2_descrip])
   ylabel('freq (Hz)')
   %xlabel('time (s)')
   subplot(4,15,6:7)
   plot(abs(mean(coh12,2)).^2,freq,'b')
   set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1])
   %set(gca,'yticklabel',[],'xticklabel',[],'ylim',rang)
   
   subplot(4,15,8:13)
   imagesc(timebase,freq,abs(pcoh_12_p_3).^2)
   axis xy
   %xlabel('time (s)')
   title([TS1_descrip,' vs ',TS2_descrip,' with ',TS3_descrip,' as predictor'])
   set(gca,'yticklabel',[],'xticklabel',[],'ylim',rang)
   subplot(4,15,14:15)
   plot(abs(mean(pcoh_12_p_3,2)).^2,freq,'b')
   set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1])
   %set(gca,'yticklabel',[],'xticklabel',[],'ylim',rang)
   
   %%%%% Signal 1 vs Signal 3
   subplot(4,15,16:20)
   imagesc(timebase,freq,abs(coh13).^2)
   axis xy
   set(gca,'ylim',rang,'xticklabel',[])
   title([TS1_descrip,' vs ',TS3_descrip])
   ylabel('freq (Hz)')
   %xlabel('time (s)')
   subplot(4,15,21:22)
   plot(abs(mean(coh13,2)).^2,freq,'b')
   set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1])
   %set(gca,'yticklabel',[],'xticklabel',[],'ylim',rang)
   
   subplot(4,15,23:28)
   imagesc(timebase,freq,abs(pcoh_13_p_2).^2)
   axis xy
   %xlabel('time (s)')
   title([TS1_descrip,' vs ',TS3_descrip,' with ',TS2_descrip,' as predictor'])
   set(gca,'yticklabel',[],'xticklabel',[],'ylim',rang)
   subplot(4,15,29:30)
   plot(abs(mean(pcoh_13_p_2,2)).^2,freq,'b')
   set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1])
   %set(gca,'yticklabel',[],'xticklabel',[],'ylim',rang)
   
   
   %%%%% Signal 2 vs Signal 3
   subplot(4,15,31:35)
   imagesc(timebase,freq,abs(coh23).^2)
   axis xy
   set(gca,'ylim',rang)
   title([TS2_descrip,' vs ',TS3_descrip])
   ylabel('freq (Hz)')
   xlabel('time (s)')
   subplot(4,15,36:37)
   plot(abs(mean(coh23,2)).^2,freq,'b')
   set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1])
   xlabel('Coherence','FontSize',7)
   %set(gca,'yticklabel',[],'ylim',rang)
   
   subplot(4,15,38:43)
   imagesc(timebase,freq,abs(pcoh_23_p_1).^2)
   axis xy
   xlabel('time (s)')
   title([TS2_descrip,' vs ',TS3_descrip,' with ',TS1_descrip,' as predictor'])
   set(gca,'yticklabel',[],'ylim',rang)
   subplot(4,15,44:45)
   plot(abs(mean(pcoh_23_p_1,2)).^2,freq,'b')
   set(gca,'yticklabel',[],'ylim',rang,'xlim',[0 1])
   xlabel(['Partial  ';'Coherence'],'FontSize',7)
   %set(gca,'yticklabel',[],'ylim',rang)
   

