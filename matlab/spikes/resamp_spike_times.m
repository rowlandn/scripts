function st_r = resamp_spike_times(cellin,resamp_rate,maxtime);

% resamp_spike_times  This functions resamples a given spike time
% array into the given frequency.  Resamp_rate is in kHz and maxtime 
% is in ms.
% 
% st_r = resampl_spike_times(spike_times,resamp_rate,maxtime);
% 
% Example: st_r = resamp_spike_times(spike_times,.4,100000);



sz=size(cellin);
len=maxtime*resamp_rate; % this number must be an integer
st_r=[];
for k=1:sz(1)
    rst=ceil(cellin{k}*resamp_rate);
    a=full(sparse(rst,1,1));
    a=resamp_rate*1e3*[a(:);zeros(len-length(a),1)];
    st_r=[st_r a];
end

    
    