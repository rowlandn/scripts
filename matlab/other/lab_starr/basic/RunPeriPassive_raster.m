function RunPeriPassive_raster(show)
% Global defines used by the other functions called by this one


global CONTIG
global VERBOSE
global RSTROWLEN

CONTIG = 50;	% Mean of 'CONTIG' adjacent points must be significant
MIN_SPK_N = 500;
RSTROWLEN = 5;  % length of each raster row in seconds
VERBOSE=false;

if ~exist('show','var')
	show = true;
end

% Edit this pattern to select time-stamp variables w/ specific names
%spkname_pattern = 'Snip\w*[abcd]';	% Kevin's pattern
spkname_pattern = '\w*';	% Accept all units
evtname_pattern = 'Accel\w*';	% Events must start w/ "Accel"


cd(uigetdir);

FileLst = dir('*.nex');
if(isempty(FileLst))
	str = pwd;
	error(['Found no NEX files in current directory - ' str ]);
end

outfid = write_text_header;

n = 0;	% Count of units processed
% For each file found in directory...
for i=1:length(FileLst)

	fname = FileLst(i).name;
	[spk_inds,spk_names] = find_nex_units( fname, spkname_pattern);
	if isempty(spk_inds)
		continue
	end

	Mvt = nex_get_events(fname,evtname_pattern);
    if isempty(Mvt)
        display(['Found no movement time stamps in file: ' fname '. Skipping.']);
        continue
    elseif length(Mvt)==1
        display(['Found only 1 mvt direction in file: ' fname '. Continuing']);
    end
	
	% For each unit in a file...
	for j=1:length(spk_inds)
		n = n+1;
		spk{n}.fname = fname;

		% save name of file & unit
		spk{n}.unitname = deblank( spk_names(j,:) );
		display(['Processing...' spk{n}.fname '....:' spk{n}.unitname]);
		
		% get spike times
        [spk_n, spk_t] = nex_ts( spk{n}.fname, spk_names(j,:),VERBOSE);
        if spk_n < MIN_SPK_N
            display(['Found only ' num2str(spk_n) ' spikes. Skipping.']);
            continue
        end
        
        % calculate isi between each mvt event
        for m = 1:Mvt(2).n
            a = Mvt(1).ts(m);
            b = Mvt(2).ts(m);
            if Mvt(1).n == Mvt(2).n && m == Mvt(2).n
                c = spk_t(end);
            else
                c = Mvt(1).ts(m+1);
            end
            avg1 = sum(spk_t >= a & spk_t <= b)/(b-a);
            avg2 = sum(spk_t >= b & spk_t <= c)/(c-b);
            write_text(outfid, fname, spk{n}.unitname, m, avg1, avg2);
        end

        
        spk{n}.raster = peripassive_raster(spk_t, RSTROWLEN);

		% for each direction of movement
		for k=1:length(Mvt)
			spk{n}.dir(k).n_reps = Mvt(k).n;
            spk{n}.dir(k).ts = Mvt(k).ts;
		end
        if( show)
            figure;
            rasterplot_passive2(spk{n},RSTROWLEN,0.9,fname(1:length(fname)-12));
        end
	end
end
fclose(outfid);

%% 
% subfunction to write text header
function outfid = write_text_header
fname = 'PeriPassive_raster.txt';
outfid = fopen(fname,'w');
fprintf(outfid,'fname\tunitname\trep_num\tMeanRateEvt1\tMeanRateEvt2\n');
return
%%
% subfunction to write a new line of text
function write_text(outfid, fname, unitname, rep_num, avg1, avg2)
fprintf(outfid, '%s\t%s\t%d\t%.3f\t%.3f\n', fname, unitname,rep_num,avg1,avg2);
return
%%
% subfunction to create raster

function raster = peripassive_raster(spk_t, RSTROWLEN)
% 
%	Make peri-event raster array of spike times - 1 row per designated row
%	length
%	NaN pads variable length rows


raster = NaN*zeros(1,1);
rst_len = size(raster,2);
rst_ncols = ceil(max(spk_t)/RSTROWLEN);
rst_start = 0;
rst_end = RSTROWLEN;
for i = 1:rst_ncols
	rstt = spk_t( find( spk_t>=rst_start & spk_t<rst_end ))-(i-1)*RSTROWLEN;
	rsttlen = size(rstt,2);
	% Resize raster array if necessary
	if rsttlen > rst_len
		raster(:,(rst_len+1:rsttlen)) = NaN;
		rst_len = size(raster,2);
	elseif rsttlen < rst_len
		rstt(rsttlen+1:rst_len) = NaN;
	end
	raster = [ raster; rstt ];
    rst_start = rst_end;
    rst_end = rst_end + RSTROWLEN;
end
% Delete 1st dummy row when finished
raster(1,:) = [];
return;

