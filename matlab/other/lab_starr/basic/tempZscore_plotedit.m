
%% Plot M1
% Display spectrogram data for relevant frequencies only (1-100Hz)
PD2plotM1 = spectrogramPD.M1Zscore(1:120,:); 
DYS2plotM1 = spectrogramDYS.M1Zscore(1:120,:);
ET2plotM1 = spectrogramET.M1Zscore(1:120,:);
 
hf1 = figure;
subplot(3,3,1);

% Lump values by std dev 
PD2plotM1(PD2plotM1 >= 2) = 2;
PD2plotM1(1<=PD2plotM1 & PD2plotM1<2) = 1;
PD2plotM1(-1< PD2plotM1 & PD2plotM1 < 1) = 0;
PD2plotM1((-2) < PD2plotM1 & PD2plotM1 <= (-1)) = (-1);
PD2plotM1(PD2plotM1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
% val1 = (min(min(PD2plotM1(:,:,:)))); 
% val2 = (max(max(PD2plotM1(:,:,:))));
% clims1 = [val1 val2];
clims1 = [-2 2];
    imagesc(taxis,faxis,PD2plotM1,clims1); 
%     colormap(gray); %makes plot black and white (greyscale)
  
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % define tick marks along y-axis
    set(gca,'YTick',(0:20:120));


hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:'); 

    % axis labels/title  
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
   
    annotation(hf1,'textbox','String',{'n=' numPD},'FontName',...
        'Arial','FontSize',9,'FitBoxToText','off','LineStyle','none',...
        'Position',[0.299 0.927 0.041 0.023]);
%     legend('boxon',['n=' num2str(numPD)],'Location',[0.299 0.927 0.041 0.023]);
        
    annotation(hf1,'textbox','String',{'Group PSD Z-scores'},'FontSize',12,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.005266 0.9639 0.1672 0.02706]);

   
    annotation(hf1,'textbox','String',{'M1'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.01874 0.8002 0.04415 0.03847]);

    annotation(hf1,'textbox','String',{'PD'},'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.211 0.948 0.058 0.043]);
%     end
    % put a color scale indicator next to the time-coherence plot
    colorbar([0.9307 0.1048 0.02354 0.8226]);

subplot(3,3,2);
% DYS2plotM1(DYS2plotM1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
% DYS2plotM1(DYS2plotM1 < (-1.96)) = (-1);
% DYS2plotM1(DYS2plotM1~=1 & DYS2plotM1~=(-1)) = 0;

DYS2plotM1(DYS2plotM1 >= 2) = 2;
DYS2plotM1(1<=DYS2plotM1 & DYS2plotM1<2) = 1;
DYS2plotM1(-1< DYS2plotM1 & DYS2plotM1 < 1) = 0;
DYS2plotM1((-2) < DYS2plotM1 & DYS2plotM1 <= (-1)) = (-1);
DYS2plotM1(DYS2plotM1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

val1 = (min(min(DYS2plotM1(:,:,:)))); 
val2 = (max(max(DYS2plotM1(:,:,:))));
% clims1 = [val1 val2];

% Plot data
imagesc(taxis,faxis,DYS2plotM1,clims1);
%     colormap(gray); %makes plot black and white (greyscale)
    
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    set(gca,'YTick',(0:20:120));
    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
   
    annotation(hf1,'textbox','String',{'n =' num2str(numDYS)},'FontName',...
        'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
        'Position',[0.579 0.927 0.041 0.023]);
%     legend('boxon',['n=' num2str(numDYS)],'Location',[0.579 0.927 0.041 0.023]);
    
    annotation(hf1,'textbox','String',{'DYS'},'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.488  0.948 0.058 0.043]);

subplot(3,3,3);
% ET2plotM1(ET2plotM1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
% ET2plotM1(ET2plotM1 < (-1.96)) = (-1);
% ET2plotM1(ET2plotM1~=1 & ET2plotM1~=(-1)) = 0;

ET2plotM1(ET2plotM1 >= 2) = 2;
ET2plotM1(1<=ET2plotM1 & ET2plotM1<2) = 1;
ET2plotM1(-1< ET2plotM1 & ET2plotM1 < 1) = 0;
ET2plotM1((-2) < ET2plotM1 & ET2plotM1 <= (-1)) = (-1);
ET2plotM1(ET2plotM1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
val1 = (min(min(ET2plotM1(:,:,:)))); % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(ET2plotM1(:,:,:))));
% clims1 = [val1 val2];

imagesc(taxis,faxis,ET2plotM1,clims1);
%     colormap(gray); %makes plot black and white (greyscale)
    
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % define Y-axis tickmarks
    set(gca,'YTick',[0:20:120]);

    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
   
    annotation(hf1,'textbox','String',{'n =' num2str(numET)},'FontName',...
        'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
        'Position',[0.859 0.927 0.041 0.023]);
%     legend('boxon',['n=' num2str(numET)],'Location',[0.859 0.927 0.041 0.023]);
    
    annotation(hf1,'textbox','String',{'ET'},'FontWeight','bold',...
    'FontSize',18,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.771  0.948 0.058 0.043]);
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
%% Plot S1
PD2plotS1 = spectrogramPD.S1Zscore(1:120,:);
DYS2plotS1 = spectrogramDYS.S1Zscore(1:120,:);
ET2plotS1 = spectrogramET.S1Zscore(1:120,:);

subplot(3,3,4);
PD2plotS1(PD2plotS1 >= 2) = 2;
PD2plotS1(1<=PD2plotS1 & PD2plotS1<2) = 1;
PD2plotS1(-1< PD2plotS1 & PD2plotS1 < 1) = 0;
PD2plotS1((-2) < PD2plotS1 & PD2plotS1 <= (-1)) = (-1);
PD2plotS1(PD2plotS1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = min(min(min(PD2plotS1(1:100,:,:))));
% val2 = max(max(max(PD2plotS1(1:100,:,:))));
% % val1 = (min(min(min(PD2plotS1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
% % val2 = (max(max(max(PD2plotS1(:,:,:)))))/10;
% clims1 = [val1 val2];
clims1 = [-2 2];

imagesc(taxis,faxis,PD2plotS1,clims1);
%     colormap(gray); %makes plot black and white (greyscale)

    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    set(gca,'YTick',(1:20:120));
    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    annotation(hf1,'textbox','String',{'S1'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.01898 0.5038 0.04415 0.03847]);

    annotation(hf1,'textbox','String',{'n =' num2str(numPDS1)},'FontName',...
        'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
        'Position',[0.299 0.628 0.041 0.023]);
%     legend('boxon',['n=' num2str(numPDS1)],'Location',[0.299 0.628 0.041 0.023]);

    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
    
subplot(3,3,5);
DYS2plotS1(DYS2plotS1 >= 2) = 2;
DYS2plotS1(1<=DYS2plotS1 & DYS2plotS1<2) = 1;
DYS2plotS1(-1< DYS2plotS1 & DYS2plotS1 < 1) = 0;
DYS2plotS1((-2) < DYS2plotS1 & DYS2plotS1 <= (-1)) = (-1);
DYS2plotS1(DYS2plotS1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = (min(min(DYS2plotS1(:,:,:)))); 
% val2 = (max(max(DYS2plotS1(:,:,:))));
% clims1 = [val1 val2];

% Plot data
imagesc(taxis,faxis,DYS2plotS1,clims1);
%     colormap(gray); %makes plot black and white (greyscale)
    
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    set(gca,'YTick',(0:20:120));
    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');

    annotation(hf1,'textbox','String',{'n =' num2str(numDYSS1)},'FontName',...
        'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
        'Position',[0.579 0.628 0.041 0.023]);
%     legend('boxon',['n=' num2str(numDYSS1)],'Location',[0.579 0.628 0.041 0.023]);

subplot(3,3,6);
ET2plotS1(ET2plotS1 >= 2) = 2;
ET2plotS1(1<=ET2plotS1 & ET2plotS1<2) = 1;
ET2plotS1(-1< ET2plotS1 & ET2plotS1 < 1) = 0;
ET2plotS1((-2) < ET2plotS1 & ET2plotS1 <= (-1)) = (-1);
ET2plotS1(ET2plotS1 <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];


% val1 = (min(min(ET2plotS1(:,:,:)))); 
% val2 = (max(max(ET2plotS1(:,:,:))));
% clims1 = [val1 val2];

imagesc(taxis,faxis,ET2plotS1,clims1);
%     colormap(gray); %makes plot black and white (greyscale)
    
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % define Y-axis tickmarks
    set(gca,'YTick',[0:20:120]);

    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    
    annotation(hf1,'textbox','String',{'n =' num2str(numETS1)},'FontName',...
        'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
        'Position',[0.859 0.628 0.041 0.023]);
%     legend('boxon',['n=' num2str(numETS1)],'Location',[0.859 0.628 0.041 0.023]);
%% Plot STN
PD2plotSTN = spectrogramPD.STNZscore(1:120,:);
DYS2plotSTN = spectrogramDYS.STNZscore(1:120,:);

subplot(3,3,7);
PD2plotSTN(PD2plotSTN >= 2) = 2;
PD2plotSTN(1<=PD2plotSTN & PD2plotSTN<2) = 1;
PD2plotSTN(-1< PD2plotSTN & PD2plotSTN < 1) = 0;
PD2plotSTN((-2) < PD2plotSTN & PD2plotSTN <= (-1)) = (-1);
PD2plotSTN(PD2plotSTN <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
% val1 = (min(min(PD2plotM1(:,:,:)))); 
% val2 = (max(max(PD2plotM1(:,:,:))));
% clims1 = [val1 val2];
clims1 = [-2 2];
    imagesc(taxis,faxis,PD2plotSTN,clims1); 
%     colormap(gray); %makes plot black and white (greyscale)
  
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % define tick marks along y-axis
    set(gca,'YTick',(0:20:120));


hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:'); 

    % axis labels/title  
    xlabel('time (sec)');
    ylabel('frequency (Hz)');

    annotation(hf1,'textbox','String',{'n =' num2str(numPDSTN)},'FontName',...
        'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
        'Position',[0.299 0.328 0.041 0.023]);
%     legend('boxon',['n=' num2str(numPDSTN)],'Location',[0.299 0.328 0.041 0.023]); 
    
    annotation(hf1,'textbox','String',{'STN'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.01898 0.1968 0.04415 0.03847]);
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
    
subplot(3,3,8);
DYS2plotSTN(DYS2plotSTN >= 2) = 2;
DYS2plotSTN(1<=DYS2plotSTN & DYS2plotSTN<2) = 1;
DYS2plotSTN(-1< DYS2plotSTN & DYS2plotSTN < 1) = 0;
DYS2plotSTN((-2) < DYS2plotSTN & DYS2plotSTN <= (-1)) = (-1);
DYS2plotSTN(DYS2plotSTN <= (-2)) = (-2);

taxis = [-2 2.48];
faxis = [0 120];

% val1 = (min(min(DYS2plotSTN(:,:,:)))); 
% val2 = (max(max(DYS2plotSTN(:,:,:))));
% clims1 = [val1 val2];

% Plot data
imagesc(taxis,faxis,DYS2plotSTN,clims1);
%     colormap(gray); %makes plot black and white (greyscale)
    
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    set(gca,'YTick',(0:20:120));
    hold;
    %plot vertical bar at movement onset
    plot([0 0],ylim,'k:');
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');

    annotation(hf1,'textbox','String',{'n =' num2str(numDYSSTN)},'FontName',...
        'Arial','FontSize',9,'FitBoxToText','on','LineStyle','none',...
        'Position',[0.579 0.328 0.041 0.023]);
%     legend('boxon',['n=' num2str(numDYSSTN)],'Location',[0.579 0.328 0.041 0.023]);
% subplot(3,3,9);
    
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);

%% Plot Coherence Z scores
Coh = false;
if Coh ~= false;
% Plot M1
PD2plotM1 = coherogramPD.M1Zscore;
DYS2plotM1 = coherogramDYS.M1Zscore;
 
hf2 = figure;
subplot(2,2,1);
% PD2plotM1(PD2plotM1 > 1.96) = 1; 
% PD2plotM1(PD2plotM1 < (-1.96)) = (-1);
% PD2plotM1(PD2plotM1~=1 & PD2plotM1~=(-1)) = 0;

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
val1 = (min(min(min(DYS2plotM1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(max(DYS2plotM1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = PD2plotM1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
   
    annotation(hf2,'textbox','String','Group coherence Z-scores',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.003469 0.97 0.1281 0.02706]);

    annotation(hf2,'textbox','String',{'M1-STN'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.003469 0.8161 0.09267 0.04213]);
    
    annotation(hf2,'textbox','String',{'PD'},'FontWeight','bold',...
    'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.2559 0.9476 0.04505 0.04335]);

%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);

subplot(2,2,2);
% DYS2plotM1(DYS2plotM1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
% DYS2plotM1(DYS2plotM1 < (-1.96)) = (-1);
% DYS2plotM1(DYS2plotM1~=1 & DYS2plotM1~=(-1)) = 0;

% val1 = min(min(min(PD2plot(1:100,:,:))));
% val2 = max(max(max(PD2plot(1:100,:,:))));
val1 = (min(min(min(DYS2plotM1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(max(DYS2plotM1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = DYS2plotM1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
   
    annotation(hf2,'textbox','String',{'DYS'},'FontWeight','bold',...
    'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.7135 0.952 0.04505 0.04335]);

%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);

% Plot S1
PD2plotS1 = coherogramPD.S1Zscore;
DYS2plotS1 = coherogramDYS.S1Zscore;

subplot(2,2,3);
% PD2plotS1(PD2plotS1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
% PD2plotS1(PD2plotS1 < (-1.96)) = (-1);
% PD2plotS1(PD2plotS1~=1 & PD2plotS1~=(-1)) = 0;
% val1 = min(min(min(PD2plotS1(1:100,:,:))));
% val2 = max(max(max(PD2plotS1(1:100,:,:))));
val1 = (min(min(min(PD2plotS1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(max(PD2plotS1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = PD2plotS1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
    annotation(hf2,'textbox','String',{'S1-STN'},'FontSize',16,...
    'FontName','Arial',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.003703 0.3454 0.09267 0.04213]);

%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
    
subplot(2,2,4);
% DYS2plotS1(DYS2plotS1 > 1.96) = 1; % temporarily limiting range of values for plotting/min/max purposes **This didn't work
% DYS2plotS1(DYS2plotS1 < (-1.96)) = (-1);
% DYS2plotS1(DYS2plotS1~=1 & DYS2plotS1~=(-1)) = 0;
% val1 = min(min(min(DYS2plotS1(1:100,:,:))));
% val2 = max(max(max(DYS2plotS1(1:100,:,:))));
val1 = (min(min(min(PD2plotS1(:,:,:)))))/10; % dividing by 10 temporarily b/c min and max values so high
val2 = (max(max(max(PD2plotS1(:,:,:)))))/10;
clims1 = [val1 val2];

    tmp1 = DYS2plotS1(1:100,:); %chopping A2plot will allow the whole colobar to be represented
    faxis_new = faxis(1:100);
    imagesc(taxis,faxis_new,tmp1,clims1);
    colormap(gray); %makes plot black and white (greyscale)
%     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
    %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:'); %6/12/09ALC: this isn't working as I adapt the script to M1 and S1 specific code, so I'm leaving it for now
%     hold(gca,'off');
    % set the y-axis direction (YDir) to have zero at the bottom
    set(gca,'YDir','normal');
    % set xlim and ylim
%     set(gca,'Xlim',[0-PRE POST]);
    set(gca,'Ylim',[0 120]);
    set (gca,'Xlim',[-2 2.5]);
    % axis labels/title
    xlabel('time (sec)');
    ylabel('frequency (Hz)');
%     annotation(hf1, 'textbox','String','PD group coherence Z-scores from M1','HorizontalAlignment','left',...
%         'Linestyle','none','Position',[ 0.003469 0.97 0.1281 0.02706]);
%     end
    % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
end
%% Below this is all old    
% %% Plot
% %faxis and taxis data are imported with time_psd data
% PD2plot = pdZscore;
% 
% % plot spectrogram for all ecog/lfp data
% hf1 = figure;
% % val1 = min(min(min(PD2plot(1:100,:,:))));
% % val2 = max(max(max(PD2plot(1:100,:,:))));
% val1 = (min(min(min(PD2plot(:,:,:))))) / 10; % dividing by 10 temporarily b/c min and max values so high
% val2 = (max(max(max(PD2plot(:,:,:))))) / 10;
% clims1 = [val1 val2];
% data_ch_names = {'e12','e23','e34','e45','e56','LFP'};
% 
% for i = 1:num_chan
%     subplot(2,3,i);
%     hold(gca,'on');
%     % make the time-frequency plot
%     tmp1 = PD2plot(1:100,:,i); %chopping A2plot will allow the whole colobar to be represented
%     faxis_new = faxis(1:100);
%     imagesc(taxis,faxis_new,tmp1,clims1);
%     colormap(gray); %makes plot black and white (greyscale)
% %     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
%     %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:');
%     hold(gca,'off');
%     % set the y-axis direction (YDir) to have zero at the bottom
%     set(gca,'YDir','normal');
%     % set xlim and ylim
% %     set(gca,'Xlim',[0-PRE POST]);
%     set(gca,'Ylim',[0 120]);
%     set (gca,'Xlim',[-2 2.5]);
%     % axis labels/title
%     xlabel('time (sec)');
%     ylabel('frequency (Hz)');
% %     if i==1
% % %         if typedata==1
% % % %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
% % %             title([outputname sprintf('\n')...
% % %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% % %                 data_ch_names{i} ' aligned to EMG Ch.' emg_ch_names{emg_i}]);
% % %         elseif typedata==2
% % % %             filename = strrep(filename,'_a.mat',''); %delete '_a.mat' from filename
% % %             title([outputname sprintf('\n')...
% % %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% % %                 data_ch_names{i} ' aligned to accel']);
% % %         else
% % %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
% %             title([outputname sprintf('\n')...
% %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% %                 data_ch_names{i} ' aligned to task button']);
% % %         end
% %     else
%         title(data_ch_names{i});
%         annotation(hf1, 'textbox','String','PD group PSD Z-scores','HorizontalAlignment','left',...
%         'Linestyle','none','Position',[ 0.003469 0.97 0.1281 0.02706]);
% %     end
%     % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
% end
%     
% %% Repeat for Dystonia data
% % Set directory path for finding relevant Dys files and import data
% %This gets directory info for each file within the DYS folder. For each file in directory, 
% %need to import Smag_mean data (eventually A2plot data) from timePSD
% %output. Then compile and average across all patients in the group.
% 
% % for DYS
% dyspath = uigetdir('', 'Select directory that contains DYS _timePSD_accel or _timePSD_emg files to be analyzed');
% dyspath = [dyspath '\'];
% cd(dyspath);
% 
% DYSdir = dir('*.mat'); % selects all .mat files in directory **may want to update to make more specific to timePSD files
% 
% numDYS = length(DYSdir);
% 
% for i = 1:numDYS
%     filename = DYSdir(i).name;
%     load(filename);
%     
%     num_chan = size(Smag_mean, 3);
% 
%     for j=1:num_chan
% % PDogram will be a large structure containing z-score analysis of
% %spectrogram and coherogram data further divided by contact 
% % PDspect will be 3D matrix of aggregated (group) spectrogram data in the format:
% %rows=freq data, col=time data, sheets = subjects
%     
%     % ------------------------------------------------------------------
%     % temporary code to normalize each subject's data to their baseline
%     % period. In the future, this will be imported along with the other
%     % variables from _timePSD_data.mat
% %     [nfchans,nframes] = size(Smag_mean(:,:,1));
% %     first = int32((0/0.05)+1); %int32(((PRE-BL(1))/t_res)+1)
% %     last = int32(1/0.05); %int32((PRE-BL(2))/t_res)
% %         for k = 1:nfchans
% %             bl = Smag_mean(k,first:last,i);
% %             blmean = mean(bl);
% %             Smag_mean(k,:,i) = Smag_mean(k,:,i)/blmean; 
% %         end
% %     %--------------------------------------------------------------------
%     DYSogram.contact_pair(j).DYSspect(:,:,i) = A2plot(:,:,j); 
%     end
% end
% 
% for i = 1:num_chan
%     % calculate mean and std deviation for spectrogram data (power) across
%     % subjects for a given contact pair
%     DYSogram.contact_pair(i).meanDYSspect = mean(DYSogram.contact_pair(i).DYSspect,3); 
%     DYSogram.contact_pair(i).stdDYSspect = std(DYSogram.contact_pair(i).DYSspect,0,3); 
%     
%     % calculate the mean of the baseline period across subjects for given
%     % contact pair
%     meanDYS_BL = mean(DYSogram.contact_pair(i).meanDYSspect(:,1:20)); % 1:20 should be replaced with variable bl from timePSD code
%     DYSogram.contact_pair(i).BLmean_value = mean(meanDYS_BL);
% %     BLgrand_mean(i,:) = BLmean_value; % accumulates mean baseline values for each contact pair
%     
%     % Want to display an error message if the baseline values are somehow
%     % not =1, but currently, this message is being displayed
%     % inappropriately (when BLmean_value = 1)
% %     if PDogram.contact_pair(i).BLmean_value ~= 1
% %        disp 'The average baseline value does not equal 1!';
% %     end
%     
%     % Calculate Z score for this contact pair
%     dysZscore(:,:,i) = (DYSogram.contact_pair(i).meanDYSspect - ...
%         DYSogram.contact_pair(i).BLmean_value) ./ DYSogram.contact_pair(i).stdDYSspect;
% 
%     zcheck = mean(dysZscore(:,1:20,i));
%     meanBLz = mean(zcheck);
% 
% end
% 
% %% Plot
% %faxis and taxis data are imported with time_psd data
% DYS2plot = dysZscore;
% 
% % plot spectrogram for all ecog/lfp data
% hf2 = figure;
% % val1 = (min(min(min(DYS2plot(1:100,:,:))))) / 10;
% % val2 = (max(max(max(DYS2plot(1:100,:,:))))) / 10;
% % clims1 = [val1 val2];
% clims1 = [-3.5 3.5]; % temporarily setting colorbar limits to these values as min/max are too big
% data_ch_names = {'e12','e23','e34','e45','e56','LFP'};
% 
% for i = 1:num_chan
%     subplot(2,3,i);
%     hold(gca,'on');
%     % make the time-frequency plot
%     tmp1 = DYS2plot(1:100,:,i); %chopping DYS2plot will allow the whole colobar to be represented
%     faxis_new = faxis(1:100);
%     imagesc(taxis,faxis_new,tmp1,clims1);
%     colormap(gray); %makes plot black and white (greyscale)
% %     imagesc(taxis,faxis,A2plot(:,:,i),clims1);
%     %plot vertical bar at movement onset
%     plot([0 0],ylim,'k:');
%     hold(gca,'off');
%     % set the y-axis direction (YDir) to have zero at the bottom
%     set(gca,'YDir','normal');
%     % set xlim and ylim
% %     set(gca,'Xlim',[0-PRE POST]);
%     set(gca,'Ylim',[0 120]);
%     set (gca,'Xlim',[-2 2.5]);
%     % axis labels/title
%     xlabel('time (sec)');
%     ylabel('frequency (Hz)');
% %     if i==1
% % %         if typedata==1
% % % %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
% % %             title([outputname sprintf('\n')...
% % %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% % %                 data_ch_names{i} ' aligned to EMG Ch.' emg_ch_names{emg_i}]);
% % %         elseif typedata==2
% % % %             filename = strrep(filename,'_a.mat',''); %delete '_a.mat' from filename
% % %             title([outputname sprintf('\n')...
% % %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% % %                 data_ch_names{i} ' aligned to accel']);
% % %         else
% % %             filename = strrep(filename,'.mat',''); %delete '.mat' from filename
% %             title([outputname sprintf('\n')...
% %                 '# epochs=' num2str(n_epochs) sprintf('\n')...
% %                 data_ch_names{i} ' aligned to task button']);
% % %         end
% %     else
%         title(data_ch_names{i});
%         annotation(hf2, 'textbox','String','DYS group PSD Z-scores','HorizontalAlignment','left',...
%         'Linestyle','none','Position',[ 0.003469 0.97 0.1281 0.02706]);
% %     end
%     % put a color scale indicator next to the time-coherence plot
%     colorbar([0.9307 0.1048 0.02354 0.8226]);
% end