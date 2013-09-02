%% m1s1 coherence program
% This program is meant to be identical to ecogCOH (previously
% recog_temp_all_graphs) to calculate coherence between two signals.
% m1s1coh calculates coherence between the designated M1 and S1 contacts
% (whereas ecogCOH estimates coherence between each ecog contact and the
% STN LFP). 

%ALC 6/10/09: The following changes need to be made to bring code up to
%date: 1. change to structure format 2. include code that limits the signal
%to 2.5 seconds. 

%% Define Constants

OFFSET = 0.5; % epoch time offset in seconds, NOTE we may change this to zero
Fs=1000;            % samples per unit time, in this case digitization at 1000 hz rate
nfft = 512; %FFT length that determines the frequencies at which the coherence is estimated. For real x and y,the length of Cxy is (nfft/2+1) if nfft is even or (nfft+1)/2 if nfft is odd. 
window = 512; %epoch data will be divided into chunks of size equal to window for coherence calculation
noverlap = 256; % how much overlap there is between windows
% nfft and window must be equal and noverlap is a fraction of those (50%,
% 75%, etc.


%% Import/Prepare file
% import filename and pathname of mat file containing ecog data created by
% apmconv7_coh
[fn pn] = uigetfile('*.mat','Select .mat containing ecog_lfp data');
cd(pn);
load([pn fn]);
% remove '_ecog_lfp.mat' ending from filename
fn = strrep(fn,'_ecog_lfp.mat','');

%% Define disease state
%Will output coh data into folders specific to PD or DYS
%Code not currently flexible enough to handle alterate diagnoses(ALC 11/25/08)
Dx = input('Enter patient diagnosis: 1=PD, 2=Dys, 3=ET, 4=Epilepsy, 5=Other : ');

% switch Dx  --attempt to make error message if wrong number entered
%     case
%         disp([num2str(input) 'is not a valid diagnosis']);
% end

%% Define ECOG contact over Motor Strip
M1_contact = input('Enter contact closest to M1 (format: [1 2 ...]): ');

%% Define ecog contact over sensory cortex
S1_contact = M1_contact -2;
if M1_contact==1
    s1 = NaN;
    menu(['There is no contact pair over M1-S1.' sprintf('\n')...
        'Function will terminate'],'OK');
    return
elseif M1_contact==2
    s1 = NaN;
    menu(['There is no contact pair over M1-S1.' sprintf('\n')...
        'Function will terminate'],'OK');
    return
end
%% Remontage Data
%such that contact 1 is referencing contact 2, 2 references 3, etc

% Recog = zeros(size(ecog_lfp_raw_data(:,end-2);
ecog1v2 = ecog_lfp_raw_data(:,1) - ecog_lfp_raw_data(:,2);
ecog2v3 = ecog_lfp_raw_data(:,2) - ecog_lfp_raw_data(:,3);
ecog3v4 = ecog_lfp_raw_data(:,3) - ecog_lfp_raw_data(:,4);
ecog4v5 = ecog_lfp_raw_data(:,4) - ecog_lfp_raw_data(:,5);
ecog5v6 = ecog_lfp_raw_data(:,5);
LFP = ecog_lfp_raw_data(:,end-2);
% %% Downsampling for Miller case (accidentally recorded at 5K)
% ecog1v2 = downsample(ecog1v2,5);
% ecog2v3 = downsample(ecog2v3,5);
% ecog3v4 = downsample(ecog3v4,5);
% ecog4v5 = downsample(ecog4v5,5);
% ecog5v6 = downsample(ecog5v6,5);
% LFP = downsample(LFP,5);

recog = [ecog1v2 ecog2v3 ecog3v4 ecog4v5 ecog5v6 LFP];


%% Round timestamps to 3 decimal places
%AlphaOmega system retains 4 decimal places, which causes rounding errors
%when parsing rest/active epochs

% as of 6/10/09 this is no longer necessary
% time_stamps = ecog_lfp_raw_data(:,end-1:end);
% time_stamps = time_stamps*1000;
% time_stamps = round(time_stamps);
% time_stamps = time_stamps/1000;
% ecog_lfp_raw_data(:,end-1:end) = time_stamps;
%% Define variables
% hf = figure;
[num_row num_col]=size(ecog_lfp_raw_data); %Need to keep raw data matrix here b/c has rest and active onsets
num_contact_pair = num_col-2; %last two columns are rest and active onsets
rest_time = ecog_lfp_raw_data(:,end-1);
active_time = ecog_lfp_raw_data(:,end);
diff = rest_time - active_time;
rest_time = rest_time(diff ~= 0);
active_time = active_time(diff~=0);
if active_time(end)== 0 %usually there will be one more rest onset than active onset, in which case we will throw out any extra zeros from active_time vector
    active_time = active_time(1:end-1);
end
num_epoch = length(active_time);

epoch_starts_r = zeros(1, num_epoch); %This will be for coherence: need to have start points of each epoch within the all_rest matrix, and will put that data here
epoch_starts_a = zeros(1,num_epoch);
epoch_ends_r = zeros(1,num_epoch);
epoch_ends_a = zeros(1,num_epoch);
%% Create empty 2D matrix where data from each epoch can be
% taken from ecog_lfp_raw_data and condensed into a single column for each
% contact pair without adding excess zeros, based on knowledge of the
% length of each epoch

epoch_length_r = zeros(1,num_epoch); %vector of length = number of epochs

%Determine length of each epoch(i) and  replace zeros in above epoch_lenghts array with length of each epoch 
for i = 1:num_epoch
        epoch_length_r(i) = int32(((active_time(i) - rest_time(i))- OFFSET) * Fs); 
end 
tot_epoch_length_r = int32(sum(epoch_length_r)); %adds all above epoch lengths together; purpose of int32 is to ensure this variable is read as an integer (has been a problem)
all_rest = zeros(tot_epoch_length_r,num_contact_pair); %creates zeros matrix: #rows= sum of all rest epochs, #columns = #contact pairs 

% Now the same for Active epochs

epoch_length_a = zeros(1,num_epoch);

for i = 1:num_epoch
        epoch_length_a(i) = int32(((rest_time(i+1) - active_time(i))- OFFSET) * Fs); 
end
tot_epoch_length_a = int32(sum(epoch_length_a));
all_active = zeros([tot_epoch_length_a num_contact_pair]);

%% Now populate those empty matrices with "raw" time domain data that is
%% segregated by activity state (rest vs movement)

for i = 1:num_contact_pair
   
    % Parse each rest/active epoch from each contact pair
    
    for j = 1:num_epoch
        
    % Determine starting point and end point of each epoch, allowing for a
    %   time offset, if needed. The start and end points must refer to a
    %   row value from the raw data matrix, thus for the first rest epoch
    %   (starts at time 0) must make alternate index in case OFFSET = 0
        start_rest(j) = int32((rest_time(j) + OFFSET) * Fs);      
        if start_rest(j) == 0
            start_rest(j) = 1; 
        end
        end_rest(j) = int32(active_time(j) * Fs - 1); %substract 1 because end of rest period is the datum just before the start of the active period that follows

    % Fill in zeros matrix with appropriate data
        if j==1
%             all_rest(1:epoch_length_r(j),i) = ecog_lfp_raw_data(start_rest(j):end_rest(j),i);
            all_rest(1:epoch_length_r(j),i) = recog(start_rest(j):end_rest(j),i);
            epoch_starts_r(j)= 1;
            epoch_ends_r(j) = epoch_length_r(j);
        else
            next_epoch_r = int32(sum(epoch_length_r(1:(j-1))) + 1); %to determine where to place real data from epochs beyond the first one, have to tell how far down column to start later epochs
            next_end_r = int32((next_epoch_r - 1) + epoch_length_r(j)); %add new epoch length to the previous epoch(s) length(s)
%             all_rest(next_epoch_r:next_end_r,i) = ecog_lfp_raw_data(start_rest(j):end_rest(j),i);
            all_rest(next_epoch_r:next_end_r,i) = recog(start_rest(j):end_rest(j),i);
            epoch_starts_r(j) = next_epoch_r; %fills epoch start times from all_rest matrix into the epoch_starts matrix for use in coherence section
            epoch_ends_r(j) = next_end_r;
        end
       
        start_active(j) = int32((active_time(j)+ OFFSET) * Fs);
        end_active(j) = int32(rest_time(j+1) * Fs - 1);
        
        if j==1
%             all_active(1:epoch_length_a(j),i) = ecog_lfp_raw_data(start_active(j):end_active(j),i);
            all_active(1:epoch_length_a(j),i) = recog(start_active(j):end_active(j),i);
            epoch_starts_a(j) = 1;
            epoch_ends_a(j) = epoch_length_a(j);
        else
            next_epoch_a = int32(sum(epoch_length_a(1:(j-1))) + 1); 
            next_end_a = int32((next_epoch_a - 1) + epoch_length_a(j));
%             all_active(next_epoch_a:next_end_a,i) = ecog_lfp_raw_data(start_active(j):end_active(j),i);
            all_active(next_epoch_a:next_end_a,i) = recog(start_active(j):end_active(j),i);
            epoch_starts_a(j) = next_epoch_a;
            epoch_ends_a(j) = next_end_a;
        end
       
    end
    
end
%% Now trying to add coherence
% this program determines the coherence of two signals, averaged using
% epochs of lengtht len_window, through the whole length of the signal

% LEN_WINDOW = 500;
% NFFT=500;
% FS=1000;

%Create empty 3D matrix where rows = signal after window averaging, 
%columns = number of comparisons (each contact vs LFP), and 
% 3rd dimesion = number of epochs
% all_epoch_meancoh_r = zeros(nfft/2+1,num_contact_pair-1,num_epoch);
% all_epoch_meancoh_a = zeros(nfft/2+1,num_contact_pair-1,num_epoch);
restm1s1coh = zeros(nfft/2+1,num_epoch);
activem1s1coh = zeros(nfft/2+1,num_epoch);
%Each epoch evaluated independently; separate graphs for each epoch

for i = 1:num_epoch
    
%     hf = figure;
    
%     Length_r(i) = epoch_length_r(i); % redundant
%     Length_a(i) = epoch_length_a(i);
%     num_loops_r(i) = floor(epoch_length_r(i)/LEN_WINDOW)-1;
%     num_loops_a(i) = floor(epoch_length_a(i)/LEN_WINDOW)-1;
%     output_r = zeros(NFFT/2+1, num_loops_r(i)); 
%     output_a = zeros(NFFT/2+1, num_loops_a(i));
    
%     First comparing all remontaged contacts to the LFP
    
%     for j = 1:num_contact_pair-1 % subtract one b/c 6 contacts means 5 contact pairs
        
        inp1_r = all_rest(epoch_starts_r(i):epoch_ends_r(i),int8(M1_contact));
        inp2_r = all_rest(epoch_starts_r(i):epoch_ends_r(i),int8(S1_contact));
        
        inp1_a = all_active(epoch_starts_a(i):epoch_ends_a(i),int8(M1_contact));
        inp2_a = all_active(epoch_starts_a(i):epoch_ends_a(i),int8(S1_contact));
        
%         for k = 1:num_loops_r(i)
%             
%             tmp1_r = inp1_r((k-1)*LEN_WINDOW+1 : k*LEN_WINDOW); 
%             tmp2_r = inp2_r((k-1)*LEN_WINDOW+1 : k*LEN_WINDOW);

% [Cxy,F] = mscohere(x,y,window,noverlap,nfft,Fs)               
            [Cxy,F] = mscohere(inp1_r,inp2_r,window,noverlap,nfft,Fs);
%             output_r(1:length(y),k) = y;
            restm1s1coh(:,i) = Cxy;
%         end
        
%         meancoh_r = mean(output_r,2);
%         all_epoch_meancoh_r(:,j,i) = meancoh_r;
%         
%         subplot(5,1,j)
%         area([8 30], [1 1],'FaceColor',[.5 .9 .6])
%         area([78 100], [1 1],'FaceColor',[.9 .85 .7])
%         plot(F, meancoh_r, 'DisplayName', 'Mean Coherence vs F,', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b');
%         hold on;
        
%         for l = 1:num_loops_a(i)
%             
%             tmp1_a = inp1_a((l-1)*LEN_WINDOW+1 : l*LEN_WINDOW);
%             tmp2_a = inp2_a((l-1)*LEN_WINDOW+1 : l*LEN_WINDOW);
            
            [Cxy,F] = mscohere(inp1_a,inp2_a,window,noverlap,nfft,Fs);
%             output_a(1:length(y),l) = y;
            activem1s1coh(:,i) = Cxy;
%         end
        
%         meancoh_a = mean(output_a,2);
%         all_epoch_meancoh_a(:,j,i) = meancoh_a;
        
%         figure(gcf);
%         h1 = subplot(5,1,j);
%         plot(F, meancoh_a, 'DisplayName', 'Mean Coherence vs F,', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Movement','color','r');
%         axis([0 150 0 1])
%         set (gca,'xtick',[0:20:150])
%         title(['Coherence of channels',num2str(j), ' and LFP in Epoch',num2str(i)]);

        
   
    
end

%% Collapse 3D matrix into 2D matrix
%This will average all the epochs together. Resulting matrix will have 
%rows = window length, columns = comparisons (ecog contact vs LFP)

mean_restm1s1coh = mean(restm1s1coh,2);
mean_activem1s1coh = mean(activem1s1coh,2);

%% Calculate inverse hyperbolic tangent for transformed coherence 
% Y = atanh(X) returns the inverse hyperbolic tangent for each element of X.

trans_coh.m1s1.rest = atanh(sqrt(mean_restm1s1coh));
trans_coh.m1s1.active = atanh(sqrt(mean_activem1s1coh));
trans_coh.freq = F;

%% Store trans coh data for contact over M1
%This will create separate fields for storing coh data from directly over
%M1. This will be handy in later code for averaging coherence between motor
%cortex and STN across subjects
trans_coh.M1 = M1_contact

%% Save and write trans coh 
%Creates .mat output file of name (filename_transcoh.mat). transcoh.mat
%files can then be used to analyze coherence data across groups
if Dx==1
    outputdir = ['C:\Users\Starr\Documents\ECOG data\Trans_coh_data\PD\'];
elseif Dx==2
    outputdir = ['C:\Users\Starr\Documents\ECOG data\Trans_coh_data\DYS\'];
elseif Dx==3
    outputdir = ['C:\Users\Starr\Documents\ECOG data\Trans_coh_data\ET\'];
elseif Dx==4
    outputdir = ['C:\Users\Starr\Documents\ECOG data\Trans_coh_data\Epilepsy\'];
else
    outputdir = ['C:\Users\Starr\Documents\ECOG data\Trans_coh_data\Other\'];
end
outputname = [outputdir fn,'_m1s1coh.mat'];
disp(['Writing ECOG/LFP  channels to:  ' outputname]);
save(outputname,'trans_coh');
% 
% % save txt file
% fname = [outputname '.txt']; %11/25/08 this is saving filename.mat.txt - can be changed later
% outfid = fopen(fname,'w');
% if(outfid == -1)
%     error(['Unable to open...' fname ]);
% end
% 
% for i = 1:num_contact_pair-1
%      fprintf(outfid, ['ECOG Ch #' num2str(i) '\t']);    
% end
% fprintf(outfid, ['LFP Ch''\t']);

%% Plotting

hf = figure('Name',(fn));
    plot(F, trans_coh.m1s1.rest, 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b','LineWidth',1.5);
    hold on;
    plot(F, trans_coh.m1s1.active, 'DisplayName', 'Active', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Movement','color','r','LineWidth',1.5);
    axis([0 150 0 1]) % wide view axis
%     axis([0 30 0 0.5]) %axis of alpha and beta freq
    set (gca,'xtick',[0:20:150])
%     daspect([20 1 1])
%     pbaspect([20 1 1])
%     p = get(h, 'pos');
%     p = p + 0.05;
%     set(h, 'pos', p);
       ylabel('Coherence')
       xlabel('F')
       title('M1-S1 Coherence');

       annotation(hf,'textbox','String',{fn},...
           'Position',[0.01135 0.9502 0.2789 0.03481],...
           'LineStyle','none', 'FitHeightToText','on');

%     set (gca,'xtick',[0:5:30])
%     set (gca,'ytick',[0:0.1:0.5])
%     title(['Coherence of channel ',num2str(i),' and LFP']);

%     hl = legend(gca,'show');
%     set(hl,'Location','SouthOutside');

% avg coherence (0 to 0.5) vs F (0 to 50)
% for i = 1:num_contact_pair-1 
%     h = subplot(4,5,(i+5));
%     plot(F, grand_meancoh_r(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b','LineWidth',1.5);
%     hold on;
%     plot(F, grand_meancoh_a(:,i), 'DisplayName', 'Active', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Movement','color','r','LineWidth',1.5);
% %     axis([0 150 0 1]) % wide view axis
%     axis([0 50 0 0.5]) %axis of alpha and beta freq
% %     set (gca,'xtick',[0:20:150])
%     set (gca,'xtick',[0:10:50])
%     set (gca,'ytick',[0:0.1:0.5])
% %     daspect([10 1 1])
%     if i==1
%        ylabel('Coherence')
%        xlabel('F')
%     end
% %     title(['Coherence of channel ',num2str(i),' and LFP']);
% end
% 
% % Transformed coherence vs F
% for i = 1:num_contact_pair-1
%     h = subplot(4,5,(i+10));
%     plot(F, trans_coh.rest(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Transformed Coherence during Rest','color','b','LineWidth',1.5);
%     hold on;
%     plot(F, trans_coh.active(:,i), 'DisplayName', 'Active', 'XDataSource', 'F', 'YDataSource', 'Transformed Coherence during Movement','color','r','LineWidth',1.5);
%     axis([0 150 0 1]) % wide view axis
% %     axis([0 30 0 0.5]) %axis of alpha and beta freq
%     set (gca,'xtick',[0:20:150])
% %     daspect([20 1 1])
%     if i==1
%        ylabel('Transformed Coherence')
%        xlabel('F')
%     end
% %     set (gca,'xtick',[0:5:50])
% %     set (gca,'ytick',[0:0.1:0.5])
% %     title(['Transformed coherence of channel ',num2str(i),' and LFP']);
% end
% 
% for i = 1:num_contact_pair-1 
%     h = subplot(4,5,(i+15));
%     plot(F, trans_coh.rest(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b','LineWidth',1.5);
%     hold on;
%     plot(F, trans_coh.active(:,i), 'DisplayName', 'Active', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Movement','color','r','LineWidth',1.5);
% %     axis([0 150 0 1]) % wide view axis
%     axis([0 50 0 1]) %axis of alpha and beta freq
% %     set (gca,'xtick',[0:20:150])
%     set (gca,'xtick',[0:10:50])
%     set (gca,'ytick',[0:0.2:1])
% %     daspect([10 1 1])
%     if i==1
%        ylabel('Transformed Coherence')
%        xlabel('F')
%     end
% %     title(['Coherence of channel ',num2str(i),' and LFP']);
% end
% %% These annotations create the contact headers at top of graph
% % Create textbox
% annotation(hf,'textbox','String',{'C1 vs LFP'},'HorizontalAlignment','center','FontWeight','bold','FitHeightToText','off','LineStyle','none',...
%     'Position',[0.1306 0.9422 0.1265 0.0485]);
% 
% % Create textbox
% annotation(hf,'textbox','HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.2816 0.935 0.1265 0.05882]);
% 
% % Create textbox
% annotation(hf,'textbox','String',{'C2 vs LFP'},...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.2867 0.9377 0.1265 0.05613]);
% 
% % Create textbox
% annotation(hf,'textbox','String',{'C3 vs LFP'},...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.4486 0.947 0.1265 0.0485]);
% 
% % Create textbox
% annotation(hf,'textbox','String',{'C4 vs LFP'},...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.6132 0.9377 0.1265 0.05613]);
% 
% % Create textbox
% annotation(hf,'textbox','String',{'C5 vs LFP'},...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FitHeightToText','off',...
%     'LineStyle','none',...
%     'Position',[0.7778 0.9501 0.1265 0.0485]);

% %% Create rectangle around contact over motor strip
% 
% % M1_contact = input('Enter contact closest to M1 (format: [1 2 ...]): ');
% 
% if M1_contact == 1
%     %if C1 over M1
%     annotation(hf,'rectangle',[0.1109 0.0462 0.1609 0.9484]);
% end
% 
% if M1_contact == 2
%     %if C2 over M1
%     annotation(hf,'rectangle',[0.2664 0.05027 0.1672 0.9484]);
% end
% 
% if M1_contact == 3
%     %if C3 over M1
%     annotation(hf,'rectangle',[0.432 0.05027 0.1672 0.9484]);
% end
% 
% if M1_contact == 4
%     %if C4 over M1
%     annotation(hf,'rectangle',[0.5914 0.04755 0.1672 0.9484]);
% end
% 
% if M1_contact ==5
%     %if C5 over M1
%     annotation(hf,'rectangle',[0.757 0.0462 0.1672 0.9484]);
% end
%% Old plotting
% for i = 1:num_contact_pair-1
%     
%     h = subplot(5,1,i);
%     plot(F, grand_meancoh_r(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b');
%     hold on;
%     plot(F, grand_meancoh_a(:,i), 'DisplayName', 'Active', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Movement','color','r');
%     axis([0 150 0 1]) % wide view axis
% %     axis([0 30 0 0.5]) %axis of alpha and beta freq
%     set (gca,'xtick',[0:20:150])
% %     set (gca,'xtick',[0:5:30])
% %     set (gca,'ytick',[0:0.1:0.5])
%     title(['Coherence of channel ',num2str(i),' and LFP']);
% end
%     hl = legend(gca,'show');
%     set(hl,'Location','SouthOutside');

% for i = 1:num_contact_pair-1
%     
%     h = subplot(5,1,i);
%     plot(F, grand_meancoh_r(:,i), 'DisplayName', 'Rest', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b');
%     hold on;
%     plot(F, grand_meancoh_a(:,i), 'DisplayName', 'Active', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Movement','color','r');
% %     axis([0 150 0 1]) % wide view axis
%     axis([0 30 0 0.5]) %axis of alpha and beta freq
% %     set (gca,'xtick',[0:20:150])
%     set (gca,'xtick',[0:5:30])
% %     set (gca,'ytick',[0:0.1:0.5])
%     title(['Coherence of channel ',num2str(i),' and LFP']);
% end
%     hl = legend(gca,'show');
%     set(hl,'Location','SouthOutside');    
    
% %% Comparing adjacent contacts (internal data check)
% for i = 1:num_epoch
%     
% %     hf = figure;
%     
%     Length_r(i) = epoch_length_r(i); % length needs to be size of epochs
%     Length_a(i) = epoch_length_a(i);
%     num_loops_r(i) = floor(Length_r(i)/LEN_WINDOW)-1;
%     num_loops_a(i) = floor(Length_a(i)/LEN_WINDOW)-1;
%     output_r = zeros(NFFT/2+1, num_loops_r(i)); 
%     output_a = zeros(NFFT/2+1, num_loops_a(i));
%     
%     for j = 1:num_contact_pair-1 - 2 
%         % num_contact_pair minus one b/c 6 contacts means 5 contact pairs,
%         % subtract an additional 2 because under remontage, we want to look
%         % at 1v2 vs 3v4, 2v3 vs 4v5, 3v4 vs 5v6
%         inp1_r = all_rest(epoch_starts_r(i):epoch_ends_r(i),j);
%         inp2_r = all_rest(epoch_starts_r(i):epoch_ends_r(i),(j+2));
%         
%         inp1_a = all_active(epoch_starts_a(i):epoch_ends_a(i),j);
%         inp2_a = all_active(epoch_starts_a(i):epoch_ends_a(i),(j+2));
%         
%         for k = 1:num_loops_r(i)
%             
%             tmp1_r = inp1_r((k-1)*LEN_WINDOW+1 : k*LEN_WINDOW); 
%             tmp2_r = inp2_r((k-1)*LEN_WINDOW+1 : k*LEN_WINDOW);
%             
%             [y,F] = mscohere(tmp1_r,tmp2_r,[],[],NFFT,FS);
%             output_r(1:length(y),k) = y;
%         end
%         
%         meancoh_r = mean(output_r,2);
%         all_epoch_meancoh_r(:,j,i) = meancoh_r;
%         
% %         subplot(5,1,j)
% %         plot(F, meancoh_r, 'DisplayName', 'Mean Coherence vs F,', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Rest','color','b');
% %         hold on;
%         
%         for l = 1:num_loops_a(i)
%             
%             tmp1_a = inp1_a((l-1)*LEN_WINDOW+1 : l*LEN_WINDOW);
%             tmp2_a = inp2_a((l-1)*LEN_WINDOW+1 : l*LEN_WINDOW);
%             
%             [y,F] = mscohere(tmp1_a,tmp2_a,[],[],NFFT,FS);
%             output_a(1:length(y),l) = y;
%         end
%         
%         meancoh_a = mean(output_a,2);
%         all_epoch_meancoh_a(:,j,i) = meancoh_a;
%         
% %         figure(gcf);
% %         subplot(5,1,j);
% %         plot(F, meancoh_a, 'DisplayName', 'Mean Coherence vs F,', 'XDataSource', 'F', 'YDataSource', 'Mean Coherence during Movement','color','r');,
% %         axis([0 150 0 1])
% %         set (gca,'xtick',[0:20:150])
% %         title(['Coherence of channels',num2str(j), ' and ', num2str(j+2),' in Epoch',num2str(i)]);
%         
%     end
%     
% end
