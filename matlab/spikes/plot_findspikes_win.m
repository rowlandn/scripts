function plot_findspikes_win(spike_times)

%spike_times = viv0416d_DCN_SC_contra_2_ant_2_lat_st;

plot(spike_times.traces); hold on
for i = 1:size(spike_times.traces,2)
    plot(spike_times.times.win{i,1}*10, spike_times.peaks{i,1},'ro')
end