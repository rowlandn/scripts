function collect_spikes = collect_spikes_SE(traces, ...
    spike_times, left, right, plotit)

% [spikes, spksavg] = collect_spikes_SE(traces, spike_times, left, right, plotit)
% 
% Gets and returns all the spikes in trace indicated by the spike_times;
% left and right are how much (in samples) you sample "around" the spike
% peaks. If plotit is 1, a figure with all spikes superimposed and averaged 
% will be displayed.
%
% <adelgado@biology.emory.edu>
%# modified by C.B.Roberts, 02/05

%spike_times = spike_times.times;

for g = 1:size(spike_times,1)
    spike_times{g,1} = spike_times{g,1}*10;
end
% assignin('base','spike_times',spike_times)

%spikes = [ ];

for h = 1:size(traces,2)
% assignin('base','h',h)
clear spikes spksavg trace
trace = traces(:,h);

    for i = 1:size(spike_times{h,1},2)
        leftplot = left;
        rightplot = right;
        if spike_times{h,1}(i)-left <= 0 
            leftplot = 0;
            rightplot = right+left;    
        end
            %# Above loop checks if the first peak is too
            %# close to the beginning.  If so, it resets 
            %# beginning ("left side") of the stored trace 
            %# to zero, and moves the right end out, so all 
            %# stored traces will be the same length.
        if spike_times{h,1}(i)+right >= length(trace) 
            trace(length(trace):length(trace)+right)=trace(end);
            
        end
            %# A similar trap to above to avoid reading over the end
            %# of the trace.  The last trace will still be a bit short,
            %# though.....
%         assignin('base','rightplot',rightplot)
%         assignin('base','leftplot',leftplot)
%         assignin('base','h',h)
        %assignin('base','i',i)
        spikes(:, i) =  trace(ceil(spike_times{h,1}(i)-leftplot): ceil(spike_times{h,1}(i)+rightplot));
        % assignin('base','spikes_collect_spikes_SE_52',spikes)
         %# The indices here are the result of some arithmetic operations,
         %# and in some special cases may not be integers.  The ceil(x)
         %# command used here rounds whatever they are to the nearest
         %# integer.  The results still seem to be "double" though....
         %# Whatever, this solves a problem of lots of warnings for 
         %# indices going from 114991 to 245628 etc.    
    end
    if size(spikes,2) == 1
        spksavg = spikes;
        collect_spk_traces{h,1} = spikes;
        collect_spk_avg{h,1} = spksavg;
    else
        spksavg = (mean(spikes'))';
        collect_spk_traces{h,1} = spikes;
        collect_spk_avg{h,1} = spksavg;
    end
end

% assignin('base','spksavg',spksavg)
% assignin('base','collect_spk_traces',collect_spk_traces)
% assignin('base','collect_spk_avg',collect_spk_avg)
        

collect_spikes.spk_traces = collect_spk_traces;
collect_spikes.spk_avg = collect_spk_avg;

if plotit == 1
    spike_x = 1:1:length(spksavg);
    spike_x = spike_x/10;
    plot(spike_x,spikes, 'k');
    hold on
    plot(spike_x,spksavg, 'g', 'LineWidth', 2);
    ylabel('Potential [mV]');
    grid on
    title('Spike Traces')
    xlabel('Time [msec]');
end
    

