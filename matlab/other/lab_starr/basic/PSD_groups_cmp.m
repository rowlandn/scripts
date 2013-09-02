

files = dir('*STN_M1_high_frq_MI.mat');

figure
hold on
x = jet(12);

Beta_Means= nan*ones(length(files),1);
Gamma_Means= nan*ones(length(files),1);
logpsd = nan*ones(length(files),257);
for i = 1 : length(files)
   
    name = files(i).name;
    load(name);
    [psdall,freq] = pwelch(ecogall(6,:),WINDOW,NOVERLAP,NFFT,Fs);
    logpsdall = log10(psdall);
    logpsd(i,:)=logpsdall;
    logpsdall([30 31 32 33 34 61 62 63 64 92 93 94 95 122 123 124 125 126 153 154 155 156 157 184 185 186 187],1)=nan;
    logmean= nanmean(logpsdall);
    plot(freq,logpsdall/logmean,'color',x(i,:))
    
    Beta_Means(i,1) = nanmean(logpsdall(8:16));
    Gamma_Means(i,1) = nanmean(logpsdall(40:52));
    
end

