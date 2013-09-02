function NS_traces = load_NS_SE(filename,traces,channels)

% load_NS_SE_03 This function loads the user-defined traces from the
% specified filename and channels into a 3-dim variable.  Output is
% {traces,trials,channels}
%
% Example: A = load_NS_SE_03('/Raw/viv06/viv0604d.data',[1:10],[1 2]);



for i = 1:size(channels,2)
    progbar = waitbar(0, 'Loading NeuroSage Traces...');
    for j = 1:size(traces,2)
        NS_traces_temp = readneurosage(filename,traces(j),channels(i));
        NS_traces(1:size(NS_traces_temp.Trials.AcqChannels.Data',1),j,i) = NS_traces_temp.Trials.AcqChannels.Data';
        waitbar(j/size(traces,2), progbar)
        disp([num2str(round((j/size(traces,2))*100)),'% complete...'])
    end
    close(progbar)
end





% assignin('base','NS_traces_temp',NS_traces_temp)
% assignin('base','i',i)

%NS_traces = NS_traces_temp.Trials.AcqChannels.Data;