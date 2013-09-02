function spikes_SE = findspikes_rms_SE1(traces,spike_times,rmst,stim_t,art_plottype)

%
% findspikes_rms_SE 
% Syntax: A = findspikes_rms_SE(traces,spike_times,pre,post,rms_thresh)
% Accepted spikes are plotted in green and rejected spikes in red 

pre = 5; %samples
post = 5; %samples

collect_spikes = collect_spikes_SE(traces,spike_times.times.win,pre,post,0);
%assignin('base','collect_spikes',collect_spikes)

for i = 1:size(traces,2)
clear rms_trace
    for j = 1:size(collect_spikes.spk_traces{i,1},2)
        err1 = collect_spikes.spk_avg{i,1}(1:pre)-collect_spikes.spk_traces{i,1}(1:pre,j);
        err2 = collect_spikes.spk_avg{i,1}(pre+1:pre+post)-collect_spikes.spk_traces{i,1}(pre+1:pre+post,j);
        err = [err1',err2']';
        sqerr = err .* err;
        rms = sqrt(mean(sqerr));
        rms_trace(j,1) = rms;
    end
    rms_all{i,1} = rms_trace;
end
%assignin('base','rms_all',rms_all)


for i = 1:size(traces,2)
clear rejected_spike_count rejected_spike_index rmst_trace
rejected_spike_count = 0;
rejected_spike_index = [];
if strmatch('auto',rmst)
    rmst_trace = 3*std(rms_all{i,1}) + mean(rms_all{i,1});
%         max_rms = max(rms_all{i,1});
%         rmst_trace = max_rms-1;
else
rmst_trace = rmst;
    end
    for j = 1:size(rms_all{i,1},1)
        if rms_all{i,1}(j) > rmst_trace
            rejected_spike_count = rejected_spike_count + 1;
            rejected_spike_index = [rejected_spike_index j];
        end       
    end
    spikes_SE.rejected_spike_counts(i,1) = rejected_spike_count;
    spikes_SE.rejected_spike_index{i,1} = rejected_spike_index;
    spikes_SE.rmst{i,1} = rmst_trace;
end
spikes_SE.rms_values = rms_all;
%assignin('base','spikes_SE',spikes_SE)

spike_count_00 = cellfun(@length,spike_times.times.win);
spike_count_01 = sum(spike_count_00);


rev_spike_times = spike_times.times.win;
%assignin('base','rev_spike_times',rev_spike_times)
rev_spike_peaks = spike_times.peaks.win;
for i = 1:size(spikes_SE.rejected_spike_index,1)
    temp_rejected_spike_times = rev_spike_times{i,1}(spikes_SE.rejected_spike_index{i,1})*10;
    find_post_stim_indices = find(temp_rejected_spike_times < stim_t*10);
    spikes_SE.rejected_spike_index{i,1}(find_post_stim_indices) = [];
    rev_spike_times{i,1}(spikes_SE.rejected_spike_index{i,1}) = [];
    rev_spike_peaks{i,1}(spikes_SE.rejected_spike_index{i,1}) = [];
end
%assignin('base','temp_rejected_spike_times',temp_rejected_spike_times)
 
spikes_SE.spike_count_01 = spike_count_01;
spikes_SE.times_win_rms = rev_spike_times;
spikes_SE.peaks_win_rms = rev_spike_peaks;
spikes_SE.collect_spikes = collect_spikes;




    

