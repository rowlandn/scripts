function burst_times = legend_salc_SE(analog_signal,spike_times,...
                                        Fs,plotit,figure_pos,sp1,sp2,sp3);

% legend_salc_SE This function identifies bursts within given spike times
% using the Poisson Surprise method of Legendy and Salcman  (J.
% Neurophysiol. 53(4):926-939, 1985) and returns the onset and offset times
% of the bursts. Spike times should be supplied as a cell containing rows
% of spike time arrays which represent individual traces (e.g., the
% output of the findspikes_win_SE function).  The function
% also returns the number of spikes (N) in each burst  (including first and
% last) and the SURPRISE value assigned to the burst. When running make
% sure you do not receive a "Log of zero" warning. This implies a burst of
% duration zero. Try supplying a spike train with a higher sampling rate
% (Fs in kHz). This function requires the functions Poisson_surprise_SE and
% FIND_TRAIN_SE to run.  If plotit is 1, the burst statistics are plotted as
% as bar graphs.  If plotit is 2, the last raw trace and corresponding burst 
% times are plotted.  If multiple graphs are to be placed on the same page, 
% use plotit = 2 and supply subplot (sp).  Else, use plotit  = 0.
%
% burst_times = legend_salc_SE(analog_signal,spike_times,Fs,plotit,sp1,sp2,sp3)
% 
% Example 1: Neuron = load_PCDX_SE('/F/Data/Raw/viv05/viv0503b.all','1-10',4);
%            spike_times = findspikes_win_SE(Neuron,10,-200,-40,.1,1,0);
%            burst_times = legend_salc_SE(Neuron,spike_times,10,1);
%             
% Example 2: Neuron = load_PCDX_SE('/F/Data/Raw/viv05/viv0503b.all','1-10',4);
%            spike_times = findspikes_win_SE(Neuron,10,-200,-40,.1,1,0);
%            burst_times = legend_salc_SE(Neuron,spike_times,10,2,2,2,1);
%                   (if you want to specify a certain subplot)

% JAG, 5/2/03  
% spike_train = A';
% spike_train = SPIKES(:,1);
% Fs = 10000;
% plotit = 1;
% Fs=50;

Fs = Fs*1000;
%spike_train=spike_times(:);
for i = 1%:size(spike_times,1)
    %assignin('base','i',i)
    spike_train=spike_times{i,1}(:);
    %assignin('base','spike_times',spike_times)
    %assignin('base','spike_train',spike_train)
    ISIs=diff(spike_train);
    %assignin('base','ISIs',ISIs)
    nu=mean(ISIs);   % mean ISI of spike train
    %assignin('base','nu',nu)
    rate=1/nu;
    %assignin('base','rate',rate)
    [ind,doublets,h]=find_train_SE(ISIs'<.5*nu,2); % finds "doublets" of short ISIs
    
k=1;
burst=[];
burst_tag=cell(1,1);
cell_cnt=0;
try % this is so it ends without an error and prints out
while k <= length(doublets)

begin_ind=doublets(k);
old_begin=begin_ind;
burst_length=2;
tent_burst=ISIs(begin_ind:(begin_ind+burst_length-1));
S_old=Poisson_surprise_SE(tent_burst,rate);
end_burst=1;
while end_burst % the burst has been initiated
   if ISIs(begin_ind+burst_length)>2*nu 
   % if the next ISI is too large then end the burst 
       end_burst=0;
       begin_ind=begin_ind+burst_length+1; 
   elseif Poisson_surprise_SE(ISIs(begin_ind:(begin_ind+burst_length)),rate)>S_old
   % otherwise if the next ISI increases the surprise then add it
       burst_length=burst_length+1;
       tent_burst=ISIs(begin_ind:(begin_ind+burst_length-1));     
   else 
   % otherwise (doesn't increase surprise and is short enough) check up to 10 intervals ahead
       counter=1;
       while counter<10
            if ISIs(begin_ind+burst_length+counter)<2*nu
            % if the next interval is small enough
               if Poisson_surprise_SE(ISIs(begin_ind:(begin_ind+burst_length+counter)),rate)>S_old 
	       % and if it increases the surprise then add all this to burst and continue    
                   burst_length=burst_length+counter+1;
                   tent_burst=ISIs(begin_ind:(begin_ind+burst_length-1)); 
                   counter=11;
               else
	       % otherwise go onto the next ISI (upto 10)
                   counter=counter+1;
               end
            else
            % if alternatively the next interval is too large then end burst
               counter=11; 
               end_burst=0;
               begin_ind=begin_ind+burst_length+1;
            end
            if counter==10
            % if gone through ten intervals and all are small but none increase the surprise
            % then end the burst
	       end_burst=0;
               begin_ind=begin_ind+burst_length+1;;
            end
       end
    end
end
back_surprise=[];
lentent=length(tent_burst);
for j=1:lentent
    back_surprise=[back_surprise Poisson_surprise_SE(tent_burst(j:end),rate)];
end
[S_max,ind]=max(back_surprise);
burst=[burst;[old_begin+ind-1 old_begin+lentent S_max]]; % record burst times and surprise value
%%assignin('base','lentent',lentent)
cell_cnt=cell_cnt+1;
burst_tag{cell_cnt}=num2str(round(100*S_max)/100);
k=min(find(doublets>=begin_ind)); % move on to next "doublet"
end
error
catch
burst_times.onset{i,1} = spike_train(burst(:,1));
burst_times.offset{i,1} = spike_train(burst(:,2));
burst_times.no_spikes{i,1} = diff(burst(:,1:2),1,2)+1;
burst_times.surp_value{i,1} = burst(:,3);
end
end
%assignin('base','burst_times',burst_times)
% Calculate mean burst frequency
for i = 1:size(burst_times.onset,1)
    total_no_bursts(i,1) = size(burst_times.onset{i,:},1);
end
%assignin('base','total_no_bursts',total_no_bursts)
mean_burst_frequency = sum(total_no_bursts)/(size(analog_signal(:),1)/Fs);
%assignin('base','Fs',Fs)
%assignin('base','analog_signal',analog_signal)
burst_times.mean_burst_frequency = mean_burst_frequency;

% Calculate mean interburst interval (ibi) cv
for j = 1:size(burst_times.onset,1)
    for k = 1:size(burst_times.onset{j,1},1)-1
            burst_intervals{j,k} = burst_times.onset{j,1}(k+1) - burst_times.offset{j,1}(k);
    end
end
for j = 1:size(burst_times.onset,1)
    std_burst_intervals{j,1} = std(cat(1,burst_intervals{j,:}));
    mean_burst_intervals{j,1} = mean(cat(1,burst_intervals{j,:}));
end
std_burst_intervals = mean(cat(1,std_burst_intervals{:}));
mean_burst_intervals = mean(cat(1,mean_burst_intervals{:}));
burst_times.ibi_cv = std_burst_intervals/mean_burst_intervals;

% Calculate mean burst duration
for j = 1:size(burst_times.onset,1)
    for k = 1:size(burst_times.onset{j,1},1)
        burst_dur{j,k} = burst_times.offset{j,1}(k) - burst_times.onset{j,1}(k);
    end
end
for j = 1:size(burst_times.onset,1)
    mean_burst_dur{j,1} = mean(cat(1,burst_dur{j,:}));
end
mean_burst_dur = mean(cat(1,mean_burst_dur{:}));
burst_times.mean_burst_duration = mean_burst_dur;

% Calculate mean intraburst spike frequency
for j = 1:size(burst_times.onset,1)
    for k = 1:size(burst_times.onset{j,1},1)
        burst_dur{j,k} = burst_times.offset{j,1}(k) - burst_times.onset{j,1}(k);
        intraburst_spike_freq{j,k} = burst_times.no_spikes{j,1}(k)/((burst_times.offset{j,1}(k) - burst_times.onset{j,1}(k))/1000);
    end
end
for j = 1:size(burst_times.onset,1)
    burst_dur_cat{j,1} = cat(1,burst_dur{j,:});
    mean_intraburst_spike_freq_cat{j,1} = mean(cat(1,intraburst_spike_freq{j,:}));
end
%assignin('base','burst_dur',burst_dur)
%assignin('base','intraburst_spike_freq',intraburst_spike_freq)
mean_intraburst_spike_freq = mean(cat(1,mean_intraburst_spike_freq_cat{:}));
burst_times.mean_intraburst_spike_freq = mean_intraburst_spike_freq;

% Calculate mean Poisson surprise
burst_times.mean_Poisson_surp = mean(cat(1,burst_times.surp_value{:}));


    
if plotit == 0
elseif plotit == 1
    
    subplot(3,6,figure_pos*2-1)
    bar(1,burst_times.mean_burst_frequency); hold on
    bar(2,burst_times.ibi_cv)
    bar(3,burst_times.mean_Poisson_surp)
    set(gca,'XTick',[1;2;3],'XTicklabel',['BF';'CV';'PS'],'FontSize',6)
    axis([0 4 0 10])
    subplot(3,6,figure_pos*2)
    bar(1,burst_times.mean_burst_duration); hold on
    bar(2,burst_times.mean_intraburst_spike_freq);
    set(gca,'XTick',[1;2],'XTicklabel',['DR';'SF'],'FontSize',6)
    axis([0 4 0 700])
    subplot(3,6,figure_pos*2-1)
    
elseif plotit == 2
    if nargin == 8
        subplot(sp1,sp2,sp3)
    end
    
    plot(analog_signal(:,end)-10,'k'); hold on
    %plot(analog_signal,'k')
    
    ind=find(burst(:,3)>2 & diff(burst(:,1:2),1,2)>2 & burst_times.offset{end,1}-burst_times.onset{end,1}<Fs);
    sp_times1=spike_train(burst(ind,1));
    sp_times2=spike_train(burst(ind,2));
    %assignin('base','spike_train',spike_train)
    plot(reshape(repmat(spike_train'*(size(analog_signal(:,end),1)/Fs),3,1),3*length(spike_train),1),repmat([0;30;NaN],length(spike_train),1)-abs(min(analog_signal(:,end)))-50,'k')
    % plot(reshape(repmat(sp_times1'*(size(analog_signal,1)/Fs),3,1),3*length(sp_times1),1),repmat([31;40;NaN],length(sp_times1),1)-abs(min(analog_signal))-20,'r')
    % plot(reshape(repmat(sp_times2'*(size(analog_signal,1)/Fs),3,1),3*length(sp_times2),1),repmat([31;40;NaN],length(sp_times2),1)-abs(min(analog_signal))-20,'k')
    %text(mean([sp_times1 sp_times2],2)'*(size(analog_signal,1)/Fs),2*ones(1,length(sp_times1))-abs(min(analog_signal))-20,burst_tag(ind),'fontsize',6)
    
%     sp_times1=spike_train(burst(:,1));
%     sp_times2=spike_train(burst(:,2));
%     
%     
%     %plot(reshape(repmat(sp_times1'*(size(analog_signal(:,end),1)/Fs),3,1),3*length(sp_times1),1),repmat([-30;0;NaN],length(sp_times1),1)-abs(min(analog_signal(:,end)))-4,'Color',[.50 .50 .50])
    plot(reshape(repmat(sp_times1'*(size(analog_signal(:,end),1)/Fs),3,1),3*length(sp_times1),1),repmat([-30;0;NaN],length(sp_times1),1)-abs(min(analog_signal(:,end)))-50,'Color',[.50 .50 .50])
    plot(reshape(repmat(sp_times2'*(size(analog_signal(:,end),1)/Fs),3,1),3*length(sp_times2),1),repmat([-30;0;NaN],length(sp_times2),1)-abs(min(analog_signal(:,end)))-50,'k')
end 
	    





