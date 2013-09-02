function [trials,neuronSpikes,neuronNames,analogRate,spikeRate] = tiProcDB(db,Nexdata,Nexvar)

events = db(3).sparse;

nTrials = size(events,2);
maxTrialLength = size(events,1);

for i=1:1:nTrials
    %if (db(1).outcome(i))
    %    fprintf('Failed trial number %i\n',i);
    %else
    %    fprintf('Good trial number %i\n',i);
    %end
    eventTimes = find(events(:,i));
    eventValues = uint32(full(events(eventTimes,i)));

    nEvents = length(eventTimes);
    
    trials(i).success = 0;

    for j=1:1:nEvents
        for k=1:1:11
            if (bitget(eventValues(j),k+1))
                switch (k)
                    case 1
                        trials(i).fix = eventTimes(j);
                    case 2
                        trials(i).hold = eventTimes(j);
                    case {3,4}
                        trials(i).target = eventTimes(j);
                        trials(i).dir = k-3;
                    case 5
                        trials(i).leave  = eventTimes(j);
                    case {6,7}
                        trials(i).touch = eventTimes(j);
                        if (trials(i).dir ~= k-6)
                            trials(i).error = 1;
                        end
                    case 8
                        trials(i).reward = eventTimes(j);
                        trials(i).success = 1;
                    case 9
                        trials(i).begin = eventTimes(j);
                    case 10
                        trials(i).end = eventTimes(j);
                end
            end
        end
    end

    trials(i).analog = squeeze(db(1).data(:,i,:));
    trials(i).spikes = squeeze(db(2).data(:,i,:));
     if (isempty(db(4).data))
        trials(i).offspikes = [];
    else
        trials(i).offspikes = squeeze(db(4).data(:,i,:));
     end
end

analogRate = db(1).AnalogSampling;
spikeRate = db(1).SpikeSampling;

if (~isempty(Nexdata) & ~isempty(Nexvar))
    neuronSpikes = Nexdata.TS(find(Nexvar.Type==0));
    neuronNames = Nexvar.Name(find(Nexvar.Type==0),:);
else
    neuronSpikes = [];
    neuronNames = [];
end
