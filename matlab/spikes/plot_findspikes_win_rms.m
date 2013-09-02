function plot_findspikes_win_rms(spike_times)

%spike_times = viv0416d_DCN_SC_contra_2_ant_2_lat_st;

plot(spike_times.traces); hold on
for i = 1:size(spike_times.traces,2)
    plot(spike_times.times.win_rms{i,1}*10, spike_times.peaks.win_rms{i,1},'go')
    plot(spike_times.times.win{i,1}(spike_times.rms_params.rejected_spike_indices{i,1})*10,spike_times.peaks.win{i,1}(spike_times.rms_params.rejected_spike_indices{i,1}),'ro')
end


% viv0416d_DCN_SC_contra_2_pos_2_lat_st.rms_params.rejected_spike_indices
% viv0416d_DCN_SC_contra_2_pos_2_lat_st.times.win{1,1}(79)