%ecogcoh determines the coherence of two signals, averaged using
%epochs of lengtht len_window, though the whole length of the signal
%
%ecogcoh_multi does the same for each of the contact pairs and displays
%them all on the same (5 subplot) graph

% 8/28/08 Not accurately plotting at least one of the graphs (#2 is clearly
% not the same as under the single runs via ecogcoh

%% 1v6
LEN_WINDOW = 500;
NFFT=500;
FS=1000;
ChA=1;
ChF=6;

Length=size(ecog_lfp_raw_data,1);
num_loops = floor(Length/LEN_WINDOW)-1;% last window omitted
output = zeros(NFFT/2+1, num_loops); 
% inp1 = ecog_lfp_raw_data(:,ChA); %selects the first input for the coherence measurement
% inp2 = ecog_lfp_raw_data(:,ChF);
inp1 = recog(:,ChA);
inp2 = recog(:,ChF);

for i = 1:num_loops
    tmp1 = inp1((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
    tmp2 = inp2((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
  
    [y,F] = mscohere(tmp1,tmp2,[],[],NFFT,FS);
    output(1:length(y),i)= y;
end
    
%NEXT PART generates mean coherence

meancoh1v6=mean(output,2); %takes the mean of each row

subplot(5,1,1)
plot(F, meancoh1v6, 'DisplayName', 'meancoh vs F', 'XDataSource', 'F', 'YDataSource', 'meancoh'); figure(gcf)
axis([0 150 0 1])
title(['coherence of channels ',num2str(ChA), ' and ', num2str(ChF)])

%% 2v6
LEN_WINDOW = 500;
NFFT=500;
FS=1000;
ChB=2;
ChF=6;

Length=size(ecog_lfp_raw_data,1);
num_loops = floor(Length/LEN_WINDOW)-1;% last window omitted
output = zeros(NFFT/2+1, num_loops); 
inp1 = ecog_lfp_raw_data(:,ChB); %selects the first input for the coherence measurement
inp2 = ecog_lfp_raw_data(:,ChF);
% inp1 = recog(:,ChB);
% inp2 = recog(:,ChF);

for i = 1:num_loops
    tmp1 = inp1((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
    tmp2 = inp2((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
  
    [y,F] = mscohere(tmp1,tmp2,[],[],NFFT,FS);
    output(1:length(y),i)= y;
end
    
%NEXT PART generates mean coherence

meancoh2v6=mean(output,2); %takes the mean of each row

subplot(5,1,2)
plot(F, meancoh, 'DisplayName', 'meancoh vs F', 'XDataSource', 'F', 'YDataSource', 'meancoh'); figure(gcf)
axis([0 150 0 1])
title(['coherence of channels ',num2str(ChB), ' and ', num2str(ChF)])

%% 3v6
LEN_WINDOW = 500;
NFFT=500;
FS=1000;
ChC=3;
ChF=6;

Length=size(ecog_lfp_raw_data,1);
num_loops = floor(Length/LEN_WINDOW)-1;% last window omitted
output = zeros(NFFT/2+1, num_loops); 
inp1 = ecog_lfp_raw_data(:,ChC); %selects the first input for the coherence measurement
inp2 = ecog_lfp_raw_data(:,ChF);
% inp1 = recog(:,ChC);
% inp2 = recog(:,ChF);

for i = 1:num_loops
    tmp1 = inp1((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
    tmp2 = inp2((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
  
    [y,F] = mscohere(tmp1,tmp2,[],[],NFFT,FS);
    output(1:length(y),i)= y;
end
    
%NEXT PART generates mean coherence

meancoh3v6=mean(output,2); %takes the mean of each row

subplot(5,1,3)
plot(F, meancoh3v6, 'DisplayName', 'meancoh vs F', 'XDataSource', 'F', 'YDataSource', 'meancoh'); figure(gcf)
axis([0 150 0 1])
title(['coherence of channels ',num2str(ChC), ' and ', num2str(ChF)])

%% 4v6
LEN_WINDOW = 500;
NFFT=500;
FS=1000;
ChD=4;
ChF=6;

Length=size(ecog_lfp_raw_data,1);
num_loops = floor(Length/LEN_WINDOW)-1;% last window omitted
output = zeros(NFFT/2+1, num_loops); 
inp1 = ecog_lfp_raw_data(:,ChD); %selects the first input for the coherence measurement
inp2 = ecog_lfp_raw_data(:,ChF);
% inp1 = recog(:,ChD);
% inp2 = recog(:,ChF);

for i = 1:num_loops
    tmp1 = inp1((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
    tmp2 = inp2((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
  
    [y,F] = mscohere(tmp1,tmp2,[],[],NFFT,FS);
    output(1:length(y),i)= y;
end
    
%NEXT PART generates mean coherence

meancoh4v6=mean(output,2); %takes the mean of each row

subplot(5,1,4)
plot(F, meancoh4v6, 'DisplayName', 'meancoh vs F', 'XDataSource', 'F', 'YDataSource', 'meancoh'); figure(gcf)
axis([0 150 0 1])
title(['coherence of channels ',num2str(ChD), ' and ', num2str(ChF)])

%% 5v6
LEN_WINDOW = 500;
NFFT=500;
FS=1000;
ChE=5;
ChF=6;

Length=size(ecog_lfp_raw_data,1);
num_loops = floor(Length/LEN_WINDOW)-1;% last window omitted
output = zeros(NFFT/2+1, num_loops); 
inp1 = ecog_lfp_raw_data(:,ChE); %selects the first input for the coherence measurement
inp2 = ecog_lfp_raw_data(:,ChF);
% inp1 = recog(:,ChE);
% inp2 = recog(:,ChF);

for i = 1:num_loops
    tmp1 = inp1((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
    tmp2 = inp2((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
  
    [y,F] = mscohere(tmp1,tmp2,[],[],NFFT,FS);
    output(1:length(y),i)= y;
end
    
%NEXT PART generates mean coherence

meancoh5v6=mean(output,2); %takes the mean of each row

subplot(5,1,5)
plot(F, meancoh5v6, 'DisplayName', 'meancoh vs F', 'XDataSource', 'F', 'YDataSource', 'meancoh'); figure(gcf)
axis([0 150 0 1])
title(['coherence of channels ',num2str(ChE), ' and ', num2str(ChF)])