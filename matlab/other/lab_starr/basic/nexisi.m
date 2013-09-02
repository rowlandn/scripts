% nexisi.m : Add-in module to be called sometime AFTER nexload.m
% requires: nexpath, nexfilename, loaded Nexfile & Nexvar data structures, 
clear neuron;

% MODIFIYABLE PARAMETERS:
if ~exist('isithreshold')
    isithreshold = 1; % minimum allowable isi, MILLISECONDS
end
if ~exist('binwidth')
    binwidth = 5;   % width of bins, MILLISECONDS
end
if ~exist('normfactor')
    normfactor = 10; % desired (normalized) mean ISI, MILLISECONDS
end

% TOGGLES USED:
if ~exist('isi_batching')
    isi_batching = 0;
end
if ~exist('isi_use_threshold')
    isi_use_threshold = 1;
end

j=0;
for i=1:Nexfile.NumVars
    if Nexvar.Type(i,:)==0
        j=j+1;
        neuron(j).file = Nexdir(n).name;
        neuron(j).name = deblank( Nexvar.Name(i,:) );
        neuron(j).ts = 1000*Nexdata.TS{i,1}; % WAS in SECONDS (in Nexfile as TICKS (*Nexfile.Frequency))
        neuron(j).isi = diff(neuron(j).ts);
        % SHOULD WE REMOVE ISI'S THAT ARE TOO SHORT ???
        if isi_use_threshold
            neuron(j).isi_threshold = isithreshold;
            clear ts_droppedindices;
            fnd = 1;
            k=0;
            while fnd  % THIS LOOP IS BASICALLY USELESS.. CUZ REALIZE WILL NEVER LOOP !
                k=k+1;
                ts_droppedindices = find(neuron(j).isi < isithreshold) +1; % ADD ONE, cuz dropping successive TS (index)
                neuron(j).ts_numdropped(k) = length(ts_droppedindices);
                if neuron(j).ts_numdropped(k) > 0
                    % Should they be appended and sorted?  i think shouldn't need >1 iteration...
                    neuron(j).ts_dropped = neuron(j).ts(ts_droppedindices); % SAVE those timestamps we are removing
                    neuron(j).ts(ts_droppedindices) = []; % REMOVE physiologically impossible timestamps from calculations
                    neuron(j).isi = diff(neuron(j).ts); % RECALCULATE ISI values.
                else
                    fnd=0; % EXIT CONDITION
                end
            end
            neuron(j).ts_numdropped = sum(neuron(j).ts_numdropped);
        end
        % ISI stats
        neuron(j).isi_count = length(neuron(j).isi);
        neuron(j).isi_min = min(neuron(j).isi);
        neuron(j).isi_max = max(neuron(j).isi);
        neuron(j).isi_mean = mean(neuron(j).isi);
        neuron(j).isi_median = median(neuron(j).isi);
        neuron(j).isi_stdev = std(neuron(j).isi);
        % Normalized ISI stats
        neuron(j).normisi = neuron(j).isi * (normfactor/neuron(j).isi_mean);
        neuron(j).normisi_min = min(neuron(j).normisi);
        neuron(j).normisi_max = max(neuron(j).normisi);
        neuron(j).normisi_mean = mean(neuron(j).normisi);
        neuron(j).normisi_median = median(neuron(j).normisi);
        neuron(j).normisi_stdev = std(neuron(j).normisi);
        % Bin ISI's
        neuron(j).isi_numbins = ceil(max(neuron(j).isi)/binwidth);
        for k=1:neuron(j).isi_numbins
            neuron(j).isi_binindices{k} = find(neuron(j).isi >= (binwidth*(k-1)) & neuron(j).isi < (binwidth*k));
            neuron(j).isi_bincounts(k) = size(neuron(j).isi_binindices{k},2); % Column dim has count
        end
        neuron(j).normisi_numbins = ceil(max(neuron(j).normisi)/binwidth);
        for k=1:neuron(j).normisi_numbins
            neuron(j).normisi_binindices{k} = find(neuron(j).normisi >= (binwidth*(k-1)) & neuron(j).normisi < (binwidth*k));
            neuron(j).normisi_bincounts(k) = size(neuron(j).normisi_binindices{k},2); % Column dim has count
        end
        % FIND THE MODE OF THE BINS: (bin with the most spikes)
        neuron(j).isi_modalvalue = max(neuron(j).isi_bincounts);
        neuron(j).isi_modalpercent = neuron(j).isi_modalvalue/neuron(j).isi_count;
        % THIS NEXT VALUE COULD BE MANY DIFFERENT BINS!  Use the FIRST, MIDDLE, or LAST one?
        bin_nums = find(neuron(j).isi_bincounts == neuron(j).isi_modalvalue);
        neuron(j).isi_modalbin = bin_nums(1); % using first now...
        neuron(j).isi_modaltime = (neuron(j).isi_modalbin-1)*binwidth+(binwidth/2);
        neuron(j).isi_burstindex = neuron(j).isi_mean/neuron(j).isi_modaltime;
        neuron(j).isi_asymmetryindex = 1/neuron(j).isi_burstindex;
        neuron(j).normisi_modalvalue = max(neuron(j).normisi_bincounts);
        neuron(j).normisi_modalpercent = neuron(j).normisi_modalvalue/neuron(j).isi_count;
        % THIS NEXT VALUE COULD BE MANY DIFFERENT BINS!  Use the FIRST, MIDDLE, or LAST one?
        bin_nums = find(neuron(j).normisi_bincounts == neuron(j).normisi_modalvalue);
        neuron(j).normisi_modalbin = bin_nums(1); % using first now...
        neuron(j).normisi_modaltime = (neuron(j).normisi_modalbin-1)*binwidth+(binwidth/2);
        neuron(j).normisi_burstindex = neuron(j).normisi_mean/neuron(j).normisi_modaltime;
        neuron(j).normisi_asymmetryindex = 1/neuron(j).normisi_burstindex;
% RLM: this is just debug info for me to use... can be commented out
disp(sprintf('################ From *.NEX file %d of %d: Neuron #%d ################', n, NUMFILES, j));
disp(sprintf('Final spike time = %g',max(neuron(j).ts)));
disp(sprintf('Smallest ISI = %g', neuron(j).isi_min));
disp(sprintf('Largest ISI = %g',neuron(j).isi_max));
neuron(j)
    end % endif type=0 
end % LOOP END: i=1:Nexfile.NumVars

if isi_batching
    % tab-delimited WRITING routine for BATCH output
    if n == 1
        fid = fopen([nexpath num2str(NUMFILES) 'Files-' basefilename(1:end-4) '_isi.txt'],'w');
        fprintf(fid,'FileName\tNeuronName\tBinWidth(ms)\t'); % Print INFO
        fprintf(fid,'ISIMean\tISIMedian\tISIStDeviation\t'); % Print ISI Stats
        fprintf(fid,'NormISIMean\tNormISIMedian\tNormISIStDeviation\t'); % Print Norm ISI Stats
        fprintf(fid,'ModalISIBin\tNormModalISIBin\t'); % Print Bin Stats
        fprintf(fid,'BurstIndex\tNormBurstIndex\t'); % Print Bin Stats
        fprintf(fid,'AsymmetryIndex\tNormAsymmetryIndex\r\n'); % Print Bin Stats
    end
    for k=1:length(neuron)
        fprintf(fid,'%s\t%s\t%g\t',Nexdir(n).name,neuron(k).name,binwidth);
        fprintf(fid,'%g\t%g\t%g\t',neuron(k).isi_mean,neuron(k).isi_median,neuron(k).isi_stdev);
        fprintf(fid,'%g\t%g\t%g\t',neuron(k).normisi_mean,neuron(k).normisi_median,neuron(k).normisi_stdev);
        fprintf(fid,'%g\t%g\t',neuron(k).isi_modaltime,neuron(k).normisi_modaltime);
        fprintf(fid,'%g\t%g\t',neuron(k).isi_burstindex,neuron(k).normisi_burstindex);
        fprintf(fid,'%g\t%g\r\n',neuron(k).isi_asymmetryindex,neuron(k).normisi_asymmetryindex);
    end
    % do not fclose(fid) unless FINAL FILE
    if n==NUMFILES
        fclose(fid);
    end
else
    fid = fopen([nexpath nexfilename(1:end-4) '_isi.txt'],'w');
    fprintf(fid,'FileName\tNeuronName\tBinWidth(ms)\t'); % Print INFO
    fprintf(fid,'ISIMean\tISIMedian\tISIStDeviation\t'); % Print ISI Stats
    fprintf(fid,'NormISIMean\tNormISIMedian\tNormISIStDeviation\t'); % Print Norm ISI Stats
    fprintf(fid,'ModalISIBin\tNormModalISIBin\t'); % Print Bin Stats
    fprintf(fid,'BurstIndex\tNormBurstIndex\t'); % Print Bin Stats
    fprintf(fid,'AsymmetryIndex\tNormAsymmetryIndex\r\n'); % Print Bin Stats
    for k=1:length(neuron)
        fprintf(fid,'%s\t%s\t%g\t',Nexdir(n).name,neuron(k).name,binwidth);
        fprintf(fid,'%g\t%g\t%g\t',neuron(k).isi_mean,neuron(k).isi_median,neuron(k).isi_stdev);
        fprintf(fid,'%g\t%g\t%g\t',neuron(k).normisi_mean,neuron(k).normisi_median,neuron(k).normisi_stdev);
        fprintf(fid,'%g\t%g\t',neuron(k).isi_modaltime,neuron(k).normisi_modaltime);
        fprintf(fid,'%g\t%g\t',neuron(k).isi_burstindex,neuron(k).normisi_burstindex);
        fprintf(fid,'%g\t%g\r\n',neuron(k).isi_asymmetryindex,neuron(k).normisi_asymmetryindex);
    end
    fclose(fid);
end
