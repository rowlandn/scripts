function RunPeriPassive(show)
% RunPeriPassive(show)
% This program tests for significant changes in firing around movement
% onset
% Created by RST, 2005-08-22
%
%	Input:
%		show - controls whether graphical output is produced for each cell
%		  default = true.
%
%	Run within a directory and this program will process all spikes in all 
%	nex files in the directory
%	


% Global defines used by the other functions called by this one
global PRE_MVT
global PST_MVT
global BIN_SZ
global CNTL_PER
global MAX_N
global ALPHA
global CONTIG

global SRCH_LO
global SRCH_HI
global VERBOSE

PRE_MVT = -1.5;
PST_MVT = 1;
BIN_SZ = 0.05;
CNTL_PER = [-1.5 -0.5];
ALPHA = 0.05;
CONTIG = 3;		% Mean of 'CONTIG' adjacent points must be significant
i=1;
j=1;
MIN_SPK_N = 500;

VERBOSE=false;

N_Pks = 2;	% max # of significant changes

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

outfid = write_text_header(N_Pks);	% Subfunction below

n = 0;	% Count of units processed
% For each file found in directory...
for i=1:length(FileLst)

	fname = FileLst(i).name;
	% Find variables of interest in file
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

		% for each direction of movement
		for k=1:length(Mvt)
			spk{n}.dir(k).n_reps = Mvt(k).n;

			[spk{n}.dir(k).histog, spk{n}.bins] = perievent_histog(spk_t, Mvt(k).ts, ...
				PRE_MVT, PST_MVT, BIN_SZ);

			spk{n}.dir(k).raster = perievent_raster(spk_t, Mvt(k).ts, ...
				PRE_MVT, PST_MVT);

			cntl_inds = find(spk{n}.bins>=CNTL_PER(1) & spk{n}.bins<CNTL_PER(2));
			test_start = max(cntl_inds)+1;
			
			[spk{n}.dir(k).chng(1),spk{n}.dir(k).cntl_mean,spk{n}.dir(k).sig_thr] = ...
				PeriEventChange(spk{n}.dir(k).histog,cntl_inds,test_start,ALPHA,CONTIG);
			
			sgn1 = spk{n}.dir(k).chng(1).sgn;
			sgn2 = [];
			if ~isempty( spk{n}.dir(k).chng(1).off_ind )
				test_start = spk{n}.dir(k).chng(1).off_ind;
				[spk{n}.dir(k).chng(2),spk{n}.dir(k).cntl_mean,q] = ...
					PeriEventChange(spk{n}.dir(k).histog,cntl_inds,test_start,ALPHA,CONTIG);
				sgn2 = spk{n}.dir(k).chng(2).sgn;
			end
			if isempty( spk{n}.dir(k).chng(1).off_ind ) | sgn2==sgn1 
				spk{n}.dir(k).chng(2).on_ind = [];
				spk{n}.dir(k).chng(2).sgn = [];
				spk{n}.dir(k).chng(2).off_ind = [];
				spk{n}.dir(k).chng(2).mean_change = [];
			end
		end
		if( show)
			make_figure(spk{n},N_Pks);
		end
				
		% write stats to file
		write_text(outfid,spk{n},N_Pks);
	end
end
fclose(outfid);
save PeriPassive spk


%------------------------------------------------------
% Subfunction to write a line of data to output file
function write_text(outfid, spk,N_Pks)
	bins = spk.bins;
	
	fprintf(outfid,'%s\t%s\t', ...
		spk.fname, spk.unitname);

	for j = 1:length(spk.dir)	% mvt directions
		mvt = spk.dir(j);
		fprintf(outfid,'%d\t%.3f\t%.3f\t',...
			 mvt.n_reps, mvt.cntl_mean, mvt.sig_thr );

		for i = 1:N_Pks
			chng = mvt.chng(i);

			% Report significant peaks
			if ~isempty(chng.on_ind)
				fprintf(outfid,'%.3f\t%.3f\t%.3f\t',...
					bins(chng.on_ind),chng.sgn,chng.mean_change);
			else
				fprintf(outfid,'-\t-\t-\t');
			end
			if ~isempty(chng.off_ind)
				fprintf(outfid,'%.3f\t',bins(chng.off_ind));
			else
				fprintf(outfid,'-\t');
			end
		end
	end
	fprintf(outfid,'\n');
	return

%------------------------------------------------------
% Subfunction to open output file and print header line
function outfid = write_text_header(N_Pks)
	fname = 'PeriPassive.txt';
	outfid = fopen(fname,'w');
	if(outfid == -1)
       error(['Unable to open...' fname ]);
	end
	fprintf(outfid,'fname\tunitname\t');
	for j = 1:2
		fprintf(outfid,'Nreps\tcntl_mean\tsig_thresh\t');
		for i = 1:N_Pks
			% For max reported signif acorr pks, freq & normalized power
			fprintf(outfid,'Onset%d\tSign%d\tMeanChange%d\tOffset%d\t',i,i,i,i);
		end
	end
	fprintf(outfid,'\n');		% EOL

	return

	

%------------------------------------------------------
% Subfunction to make figure of results
function make_figure(s, N_Pks)
	bins = s.bins;
	
	%%%%%%%%%%%%%% Plotting
	% Set up axes
	MARGIN = 0.06;	
	TOP = 1-MARGIN;		% Top margin of page
	LEFT = MARGIN;	% Left margin of page
	WIDTH = (1-3*MARGIN)/2;	% give space for 3 margin widths including middle
	HEIGHT = 0.4;	% Height of histograms
	figure
	set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
	% Size to make it look good
	c = get(gcf);
	c.Position(2) = 275;
	c.Position(3) = 870;
	c.Position(4) = 680;
	set(gcf,'Position',c.Position);

	% Find max across directions
	for j = 1:length(s.dir)
		x(j) = max(s.dir(j).histog);
		ymax = max(x)+5;
	end

	% Plot for 2 movements
	for j=1:length(s.dir)
		mvt = s.dir(j);
		left = MARGIN + (WIDTH+MARGIN)*(j-1);
		width = WIDTH;
		height = HEIGHT;      
		bottom = TOP-height;
		subplot('position',[left bottom width height]);

		h=bar(bins, mvt.histog,'histc');
		set(h,'FaceColor',[0.4,0.4,0.4],'EdgeColor','k');
		xlim([min(bins) max(bins)]);
		ylim([0 ymax]);
		ylm = ylim;
		xlm = xlim;
		hold on
		plot(xlim,[mvt.cntl_mean mvt.cntl_mean],'k-');
		plot(xlim,[mvt.cntl_mean+mvt.sig_thr mvt.cntl_mean+mvt.sig_thr],'k:');
		plot(xlim,[mvt.cntl_mean-mvt.sig_thr mvt.cntl_mean-mvt.sig_thr],'k:');
		xlabel('seconds');
		ylabel('spikes/sec');

		if isempty(mvt.chng(1).on_ind)
			text( mean([max(bins) min(bins)]), ylm(2)/2, 'No sig change found',...
					'HorizontalAlignment','center','Color','r');
		else
			for i = 1:2
				chng = mvt.chng(i);

				if ~isempty(chng.on_ind);
					plot([bins(chng.on_ind) bins(chng.on_ind)],ylm,'r:')
				end
				if ~isempty(chng.off_ind);
					plot([bins(chng.off_ind) bins(chng.off_ind)],ylm,'r:')
				end
			end
		end

		
		if j ==1
			title([ s.fname ':    ' s.unitname],'Interpreter','none','FontSize',14,...
				'Position',[xlm(2),ylm(2)+5,0]);
		end

		bottom = bottom-height-MARGIN;
		subplot('position',[left bottom width height]);
		rasterplot(mvt.raster,0.9,'');
	end

return


