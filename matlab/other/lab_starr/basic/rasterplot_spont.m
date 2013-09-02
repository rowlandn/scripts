function rasterplot_spont(show)
% Global defines used by the other functions called by this one


global CONTIG
global VERBOSE
global RSTROWLEN

CONTIG = 50;	% Mean of 'CONTIG' adjacent points must be significant
RSTROWLEN = 5;  % length of each raster row in seconds
VERBOSE=false;

if ~exist('show','var')
	show = true;
end

% Edit this pattern to select time-stamp variables w/ specific names
spkname_pattern = '\w*';	% Accept all units


cd(uigetdir);

FileLst = dir('*.nex');
if(isempty(FileLst))
	str = pwd;
	error(['Found no NEX files in current directory - ' str ]);
end


n = 0;	% Count of units processed
% For each file found in directory...
for i=1:length(FileLst)

	fname = FileLst(i).name;
	[spk_inds,spk_names] = find_nex_units( fname, spkname_pattern);
	if isempty(spk_inds)
		continue
    end
	
	% For each unit in a file...
	for j=1:length(spk_inds)
		n = n+1;
		spk{n}.fname = fname;

		% save name of file & unit
		spk{n}.unitname = deblank( spk_names(j,:) );
		% display(['Processing...' spk{n}.fname '....:' spk{n}.unitname]);
		
		% get spike times
        [spk_n, spk_t] = nex_ts( spk{n}.fname, spk_names(j,:),VERBOSE);

        
        % subfunction below
        spk{n}.raster = raster_spont(spk_t, RSTROWLEN);

        if( show)
            raster_spont_figure(spk{n},RSTROWLEN,0.9,fname(1:length(fname)-12))
        end
	end
end

%%
% subfunction to create raster

function raster = raster_spont(spk_t, RSTROWLEN)
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
%%
function raster_spont_figure(s,RSTROWLEN,tickheight,txt,xlm,color,marker)

figure;

spike_mat = s.raster;

if ~exist('xlm','var')
	xlm = [min(min(spike_mat)) max(max(spike_mat))];
end
if ~exist('color','var')
    f = get(gcf);
    if sum(f.Color) == 0
        color = 'w';
    else
        color = 'k';
    end
end

if ~exist('marker','var')
    marker.dir(1) = 'w';
    marker.dir(2) = 'k';
end


missing = find( isnan(spike_mat(:,1)) );
spike_mat(missing,:) = [];

[nrows,maxspks] = size(spike_mat);


trial_n = zeros(nrows,maxspks);

for i = 1:nrows
	trial_n(i,:) = i;
end

% transform spike times & trial numbers into row vectors
x = reshape(spike_mat',1,[]);
y = reshape(trial_n',1,[]);

newx = zeros(2,length(x));
newx(1,:) = x;
newx(2,:) = x;
newy = zeros(2,length(y));
newy(1,:) = y;
newy(2,:) = y + tickheight;

h = plot(newx,newy, 'LineWidth', 0.01,'color', color);

ylim([tickheight nrows+tickheight]); 
x = get(gca);
ind = find( rem(x.YTick,1)==0 );
ytk = x.YTick(ind);
set(gca,'YTick',ytk);

xlim( xlm );
title(['Unitname: ' txt]);
