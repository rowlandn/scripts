function dataplot
% This function plots the DBS OFF 30 sec rest data with the DBS ON rest
% data. Input comes from the output of the ecogPSDCOHrest function. Puts
% out a .fig file as well as a .mat file with modified 'psdall' variables. 

%% initialize variables
XLIM_SPEC1 = [0 150];     % range of spectral frequency for plotting (zoomed out)
XLIM_SPEC2 = [0 50];        % range of spectral frequency for plotting (zoomed in)
YLIM_NORM = [0 1.3];    % ylim for plotting normalized PSD
YLIM_LOG = [-3 3];      % ylim for for plotting log PSD
%---frequency ranges for graphing---
FREQ_LO = [8 30];       % beta band
% FREQ_MED = [35 57];     %low gamma band
FREQ_HI = [78 100]; 

FREQ_QPSD = [4 13;...     % delta alpha band
             13 22;...   % low beta band
             22 31;...   % high beta band
             31 55;...   % low gamma band
             76 100];    % high gamma band  
%% Get files
[fn pn] = uigetfile('*_ecogPSDrest.mat','Select baseline file');
cd(pn);
% load ecog data
load([pn fn]);
filename1 = fn(1:6);
DBS = [onoff];

        basepsdall = psdall;
        if size(basepsdall,2) > 5
        basepsdall(:,6) = []; % if we're using a 30 sec rest recording here, we don't
        %need the LFP data
        end
        basenormpsdall = zeros(size(basepsdall)); %initialize normpsdall
    
[fn pn] = uigetfile('*_ecogPSDrest.mat','Select comparison file');
cd(pn);
% load ecog data
load([pn fn]);
filename2 = fn(1:6);
DBS = [DBS onoff];         
        
        varpsdall = psdall; 
        if size(varpsdall,2) > 5
        varpsdall(:,6) = [];
        end
        varnormpsdall = zeros(size(varpsdall));
   

        %norm_idx=find(freq>8 & freq<100); % use norm_idx to normalize by max power between 8-100Hz, SAS 11/24/09
        betanorm_idx = find(freq>8 & freq<35); 
for i=1:size(basepsdall,2)
        %     normpsdall(:,i)=psdall(:,i)/max(psdall(:,i)); % normalize each column to its max value
        basenormpsdall(:,i)=basepsdall(:,i)/max(basepsdall(betanorm_idx(1):betanorm_idx(end),i)); % normalize each column to its max value
        varnormpsdall(:,i)=varpsdall(:,i)/max(basepsdall(betanorm_idx(1):betanorm_idx(end),i));
end
        % log PSD
        baselogpsdall = log10(basepsdall);
        varlogpsdall = log10(varpsdall);
        diff = varlogpsdall - baselogpsdall;
        
 outputname = [filename1 'vs' filename2];  
 save([outputname '_ecogPSDrest.mat'], 'basepsdall', 'varpsdall', 'basenormpsdall', 'varnormpsdall',...
     'psdall', 'baselogpsdall', 'varlogpsdall', 'diff');
     
        

%% Plot 
hf1 = figure;
for i = 1:5
    % 1st row subplots
    subplot(4,5,i);      
    ha = gca;
    hold(ha,'on');
    plot(freq, basenormpsdall(:,i),...
        freq, varnormpsdall(:,i),'r','LineWidth',2);              
        title(['e' num2str(i) num2str(i+1)]);        
        if i == 1            
            xlabel('frequency (Hz)');
            ylabel('normalized PSD');           
             if DBS == ['nn'] 
            h2 = legend('Base','Comp');
             else
            h2 = legend('DBS OFF','DBS ON');
             end        
            set(h2,'FontSize',5,'Box','off');            
        end 
       
    set(ha,'YLim',YLIM_NORM);
    set(ha,'XLim',XLIM_SPEC1);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');
    
    %second row subplots
    subplot(4,5,5+i);
    ha = gca;       % create handle for current axes
    hold(ha,'on');
    plot(freq, basenormpsdall(:,i),...
        freq, varnormpsdall(:,i),'r','LineWidth',2);
    if i == 1            
            xlabel('frequency (Hz)');
            ylabel('normalized PSD');           
             if DBS == ['nn'] 
            h2 = legend('Base','Comp');
             else
            h2 = legend('DBS OFF','DBS ON');
             end        
            set(h2,'FontSize',5,'Box','off');            
        end 
    set(ha,'YLim',YLIM_NORM);
    set(ha,'XLim',XLIM_SPEC2);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma
    hold(ha,'off');
    
     % 3rd row subplots
    subplot(4,5,2*5+i);
    ha = gca;
    hold(ha,'on');
    plot(freq,baselogpsdall(:,i),...
        freq,varlogpsdall(:,i),...
        'r','LineWidth',2);
    if i == 1            
            xlabel('frequency (Hz)');
            ylabel('normalized PSD');           
             if DBS == ['nn'] 
            h2 = legend('Base','Comp');
             else
            h2 = legend('DBS OFF','DBS ON');
             end        
            set(h2,'FontSize',5,'Box','off');            
        end 
    set(ha,'YLim',YLIM_LOG);
    set(ha,'XLim',XLIM_SPEC1);
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); 
    hold(ha,'off');
    
     % 4th row subplot
    subplot(4,5,3*5+i);
    ha = gca;
    hold(ha,'on');
    plot(freq, diff(:,i),'r','LineWidth',2);
    plot(XLIM_SPEC1,[0 0],'--k');
    if i == 1
        ylabel(['difference of' sprintf('\n')...
            'log PSD' sprintf('\n')...
            'between ON and OFF'],...
            'FontSize',8);
        xlabel('frequency (Hz)');
        hl = legend('ON');
        set(hl,'FontSize',5,'Box','off');
    end
    set(ha,'XLim',XLIM_SPEC1);
    set(ha,'YLim',[-2 1]); 
    FFT_band_fill(ha,FREQ_HI,FREQ_LO); % fill in colors for beta and gamma bands
    hold(ha,'off');

end


saveas(hf1,[outputname  '_ecogPSDrest'],'fig');
    
end
    
  



function FFT_band_fill(ha,FREQ_HI,FREQ_LO)
% FFT_band_fill fills bands of specified frequncies on axes specified by
% the handle ha

ylm = get(ha,'YLim');
x_lo = [FREQ_LO(1) FREQ_LO(2) FREQ_LO(2) FREQ_LO(1)];
x_hi = [FREQ_HI(1) FREQ_HI(2) FREQ_HI(2) FREQ_HI(1)];
y = [ylm(1) ylm(1) ylm(2) ylm(2)];
fill(x_lo,y,'g',...
    'EdgeColor','none',...
    'FaceAlpha',0.3);
fill(x_hi,y,'y',...
    'EdgeColor','none',...
    'FaceAlpha',0.3);
return;
end