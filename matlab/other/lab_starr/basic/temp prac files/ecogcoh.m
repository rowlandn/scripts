%this program determines the coherence of two signals, averaged using
%epochs of lengtht len_window, though the whole length of the signal

LEN_WINDOW = 500;
NFFT=500;
FS=1000;
ChA=1;
ChB=6;

% try taking ecog_lfp_raw_data and remontaging to make a matrix where the
% columns represent: column 1 bipolar ecog 1-2, column 2 bipolar ecog 2-3
% try having the program plot all 6 coherence measurements on one plot,
% like sho's ecog fft

Length=size(ecog_lfp_raw_data,1);
num_loops = floor(Length/LEN_WINDOW)-1;% last window omitted
output = zeros(NFFT/2+1, num_loops); 
% inp1 = ecog_lfp_raw_data(:,ChA); %selects the first input for the coherence measurement
% inp2 = ecog_lfp_raw_data(:,ChB);
% inp1 = recog(:,ChA); % recog=matrix of remontaged ecog data
% inp2 = recog(:,ChB);
inp1 = allrest(:,ChA); %allrest=matrix of ecog data during rest blocks
inp2 = allrest(:,ChB);
% inp1 = allactive(:,ChA);
% inp2 = allactive(:,ChB);

for i = 1:num_loops
    tmp1 = inp1((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
    tmp2 = inp2((i-1)*LEN_WINDOW+1 : i*LEN_WINDOW);
  
    [y,F] = mscohere(tmp1,tmp2,[],[],NFFT,FS);
    output(1:length(y),i)= y;
end
    
%NEXT PART generates mean coherence

meancoh=mean(output,2); %takes the mean of each row

plot(F, meancoh, 'DisplayName', 'meancoh vs F', 'XDataSource', 'F', 'YDataSource', 'meancoh'); figure(gcf)
axis([0 150 0 1])
set (gca,'xtick',[0:20:150])
title(['coherence of channels ',num2str(ChA), ' and ', num2str(ChB)])

