function delta = delta_SE(spike_times,max_time,Fs)

% DELTA_SE This function convolves spike times with a  
% delta  function, which is simply 1/dt. Enter the 
% length of the desired  trace (max_time) in ms. Enter 
% sampling frequency (Fs) in kHz.
%
% delta = delta_SE(spike_times,max_time,Fs)
%
% Example 1: A = load_PCDX_SE('/Raw/viv05/viv0518b.all','1-10',4);
%            spike_times = findspikes_win_SE(A,10,{-300 -50 .1 2},1);
%            delta = delta_SE(spike_times,10000,10);

% Scale parameters
max_data = Fs*max_time;
dt = 1/(Fs*1000); 

if size(spike_times,2) == 2
    spike_times = spike_times(:,1);
else
end

delta = zeros(max_data,size(spike_times,1));

for i = 1:size(spike_times,1)
    spike_time_index = round(spike_times{i,1}*Fs);
    delta(spike_time_index,i) = 1/dt;
end