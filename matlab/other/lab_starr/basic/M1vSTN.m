% M1vSTN plots data from M1 contact vs STN lfp

%% Part I
% Plot allfreq box 1 (freq with max log power) for all PD patients ("ARM)
% as "proof of code"

% grab PD data
pdpath = uigetdir('', 'Select directory that contains PD _ecogPSD files to be analyzed');
pdpath = [pdpath '\'];
cd(pdpath);

PDdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
numPD = length(PDdir);

PDallrest = [];
counterPD = 0;

for i=1:numPD
    filename = PDdir(i).name;
    load(filename);
    PDallrest(i,1) = allfreq(1,1,1);%#ok<AGROW> 
    PDallrest(i,2) = allfreq(1,1,3);%#ok<AGROW> 
    PDallactive(i,1) = allfreq(2,1,1);%#ok<AGROW> 
    PDallactive(i,2) = allfreq(2,1,3);%#ok<AGROW> 
    counterPD = counterPD+1;
end

[R,P] = corrcoef(PDallrest);
R = R(2,1);
P = P(2,1);

[Ra Pa] = corrcoef(PDallactive);
Ra = Ra(2,1);
Pa = Pa(2,1);

h=figure;
h=subplot(2,1,1);
scatter(PDallrest(:,1), PDallrest(:,2),'MarkerEdgeColor','b')
axis([0 30 0 30]);
title('PD allfreq ARM');
xlabel({'M1'});
ylabel({'STN'});
legend(['r^2 = ' num2str(R) '; p = ' num2str(P)], 'Location', 'EastOutside'); 

h=subplot(2,1,2);
scatter(PDallactive(:,1), PDallactive(:,2),'MarkerEdgeColor','r')
axis([0 30 0 30]);
legend(['R = ' num2str(Ra) '; p = ' num2str(Pa)], 'Location', 'EastOutside'); 
xlabel({'M1'});
ylabel({'STN'});

%% Part II
% since Part I works: move on to logrest
%Note, this code works, but logrest does not produce a correlation between
%m1 and stn
% normalizing produces some correlation
PDrNORMm1= [];
PDrNORMSTN = [];
for i=1:numPD
    filename = PDdir(i).name;
    load(filename);
    
    restm1 = rest(order(1,1),4:53);
    NORMrestM1 = restm1./sum(restm1);
    restSTN = rest(order(1,3),4:53);
    NORMrestSTN = restSTN./sum(restSTN);
    
    PDrNORMm1 = [PDrNORMm1  NORMrestM1]; 
    PDrNORMSTN = [PDrNORMSTN NORMrestSTN]; 
    counterPD = counterPD+1;
end
PDrNORMm1 = PDrNORMm1';
PDrNORMSTN = PDrNORMSTN';
PDrest = [PDrNORMm1 PDrNORMSTN];

[R,P] = corrcoef(PDrest);
R = R(2,1);
P = P(2,1);

H2 = figure('Name','M1vsSTN in PD, rest');
scatter(PDrest(:,1), PDrest(:,2))
legend(['r^2 = ' num2str(R) '; p = ' num2str(P)]); 
% hold on;
% xplot1 = linspace(0, 0.35); %linspace(xlim(1), xlim(2))
% lin_reg = polyfit(PDrNORMm1, PDrNORMSTN,1);
% yplot1 = polyval(lin_reg, xplot1);
% plot(xplot1,yplot1,'color','black');
% title('M1vsSTN in PD, rest');
%%
% repeats above code, but now for active condition
PDaNORMm1= [];
PDaNORMSTN = [];
for i=1:numPD
    filename = PDdir(i).name;
    load(filename);
    
    activem1 = active(order(1,1),4:53);
    NORMactiveM1 = activem1./sum(activem1);
    activeSTN = active(order(1,3),4:53);
    NORMactiveSTN = activeSTN./sum(activeSTN);
    
    PDaNORMm1 = [PDaNORMm1  NORMactiveM1]; 
    PDaNORMSTN = [PDaNORMSTN NORMactiveSTN]; 
    counterPD = counterPD+1;
end
PDaNORMm1 = PDaNORMm1';
PDaNORMSTN = PDaNORMSTN';
PDactive = [PDaNORMm1 PDaNORMSTN];

[R,P] = corrcoef(PDactive);
Ra = R(2,1);
Pa = P(2,1);

H3 = figure('Name','M1vsSTN in PD, active');
scatter(PDactive(:,1), PDactive(:,2),'MarkerEdgeColor','r')
legend(['r^2 = ' num2str(Ra) '; p = ' num2str(Pa)]); 
hold on;
xplot1 = linspace(0, 0.35); %linspace(xlim(1), xlim(2))
lin_reg = polyfit(PDrNORMm1, PDrNORMSTN,1);
yplot1 = polyval(lin_reg, xplot1);
plot(xplot1,yplot1,'color','black');
title('M1vsSTN in PD, active');

%% Plotting M1 vs STN in Beta band
%Specifically, plotting the %Power in Beta in M1 vs the 
%Power in Beta in STN

PDbetarest = [];
counterPD = 0;

for i=1:numPD
    filename = PDdir(i).name;
    load(filename);
%     if isnan(order(1,3));
%         continue;
%     end 
    PDbetarest(i,1) = (subfreq(2,2,1)+subfreq(3,2,1));%#ok<AGROW> 
    PDbetarest(i,2) = (subfreq(2,2,3)+subfreq(3,2,3));%#ok<AGROW> 
    PDbetaactive(i,1) = (subfreq(2,4,1)+subfreq(3,4,1));%#ok<AGROW> 
    PDbetaactive(i,2) = (subfreq(2,4,3)+subfreq(3,4,3));%#ok<AGROW>
    counterPD = counterPD+1;
end

PDbetarest=100*PDbetarest;
PDbetaactive=100*PDbetaactive;

[R,P] = corrcoef(PDbetarest);
R = R(2,1);
P = P(2,1);

[Ra Pa] = corrcoef(PDbetaactive);
Ra = Ra(2,1);
Pa = Pa(2,1);

h4=figure;
h2=subplot(2,1,1);
scatter(PDbetarest(:,1), PDbetarest(:,2),'MarkerEdgeColor','b')
axis([0 100 0 100]);
title('PD Beta %Power ARM');
xlabel({'M1'});
ylabel({'STN'});
legend(['r^2 = ' num2str(R) '; p = ' num2str(P)], 'Location', 'EastOutside'); 
% hold on;
% xplot1 = linspace(0, 100); %linspace(xlim(1), xlim(2))
% lin_reg = polyfit(PDrNORMm1, PDrNORMSTN,1);
% yplot1 = polyval(lin_reg, xplot1);
% plot(xplot1,yplot1,'color','black');

h2=subplot(2,1,2);
scatter(PDbetaactive(:,1), PDbetaactive(:,2),'MarkerEdgeColor','r')
axis([0 100 0 100]);
legend(['r^2 = ' num2str(Ra) '; p = ' num2str(Pa)], 'Location', 'EastOutside'); 
xlabel({'M1'});
ylabel({'STN'});
hold on;
xplot1 = linspace(0, 100); %linspace(xlim(1), xlim(2))
lin_reg = polyfit(PDrNORMm1, PDrNORMSTN,1);
yplot1 = polyval(lin_reg, xplot1);
plot(xplot1,yplot1,'color','black');
%% Part III - dystonia cases
% grab DYS data
dyspath = uigetdir('', 'Select directory that contains DYS _ecogPSD files to be analyzed');
dyspath = [dyspath '\'];
cd(dyspath);

DYSdir = dir('*_ecogPSD.mat'); % selects '*_ecogPSD.mat files' that were outputs of the ecogPSD4quant analysis script
numDYS = length(DYSdir);

DYSallrest = [];
counterDYS = 0;

for i=1:numDYS
    filename = DYSdir(i).name;
    load(filename);
    DYSallrest(i,1) = allfreq(1,1,1);%#ok<AGROW> 
    DYSallrest(i,2) = allfreq(1,1,3);%#ok<AGROW> 
    DYSallactive(i,1) = allfreq(2,1,1);%#ok<AGROW> 
    DYSallactive(i,2) = allfreq(2,1,3);%#ok<AGROW>
    counterDYS = counterDYS+1;
end

[R,P] = corrcoef(DYSallrest);
R = R(2,1);
P = P(2,1);

[Ra Pa] = corrcoef(DYSallactive);
Ra = Ra(2,1);
Pa = Pa(2,1);

h4=figure;
h2=subplot(2,1,1);
scatter(DYSallrest(:,1), DYSallrest(:,2),'MarkerEdgeColor','b')
axis([0 30 0 30]);
title('DYS allfreq ARM');
xlabel({'M1'});
ylabel({'STN'});
legend(['R = ' num2str(R) '; p = ' num2str(P)], 'Location', 'EastOutside'); 

h2=subplot(2,1,2);
scatter(DYSallactive(:,1), DYSallactive(:,2),'MarkerEdgeColor','r')
axis([0 30 0 30]);
legend(['R = ' num2str(Ra) '; p = ' num2str(Pa)], 'Location', 'EastOutside'); 
xlabel({'M1'});
ylabel({'STN'}); 
%% Using PSD data (freq 1-100Hz)
DYSrNORMm1= [];
DYSrNORMSTN = [];
for i=1:numDYS
    filename = DYSdir(i).name;
    load(filename);
    
    restm1 = rest(order(1,1),4:53);
    NORMrestM1 = restm1./sum(restm1);
    if isnan(order(1,3));
        continue;
    end 
    restSTN = rest(order(1,3),4:53);
    NORMrestSTN = restSTN./sum(restSTN);
    
    DYSrNORMm1 = [DYSrNORMm1  NORMrestM1]; 
    DYSrNORMSTN = [DYSrNORMSTN NORMrestSTN]; 
    counterDYS = counterDYS+1;
end
DYSrNORMm1 = DYSrNORMm1';
DYSrNORMSTN = DYSrNORMSTN';
DYSrest = [DYSrNORMm1 DYSrNORMSTN];

[R,P] = corrcoef(DYSrest);
R = R(2,1);
P = P(2,1);

H5 = figure('Name','M1vsSTN in DYS, rest');
scatter(DYSrest(:,1), DYSrest(:,2))
legend(['r^2 = ' num2str(R) '; p = ' num2str(P)]); 
% hold on;
% xplot1 = linspace(0, 0.35); %linspace(xlim(1), xlim(2))
% lin_reg = polyfit(PDrNORMm1, PDrNORMSTN,1);
% yplot1 = polyval(lin_reg, xplot1);
% plot(xplot1,yplot1,'color','black');
% title('M1vsSTN in DYS, rest');
% %%
%% repeats above code, but now for active condition
DYSaNORMm1= [];
DYSaNORMSTN = [];
for i=1:numDYS
    filename = DYSdir(i).name;
    load(filename);
    
    activem1 = active(order(1,1),4:53);
    NORMactiveM1 = activem1./sum(activem1);
    if isnan(order(1,3));
        continue;
    end
    activeSTN = active(order(1,3),4:53);
    NORMactiveSTN = activeSTN./sum(activeSTN);
    
    DYSaNORMm1 = [DYSaNORMm1  NORMactiveM1]; 
    DYSaNORMSTN = [DYSaNORMSTN NORMactiveSTN]; 
    counterDYS = counterDYS+1;
end
DYSaNORMm1 = DYSaNORMm1';
DYSaNORMSTN = DYSaNORMSTN';
DYSactive = [DYSaNORMm1 DYSaNORMSTN];

[R,P] = corrcoef(DYSactive);
Ra = R(2,1);
Pa = P(2,1);

H6 = figure('Name','M1vsSTN in DYS, active');
scatter(DYSactive(:,1), DYSactive(:,2),'MarkerEdgeColor','r')
legend(['r^2 = ' num2str(Ra) '; p = ' num2str(Pa)]); 
% hold on;
% xplot1 = linspace(0, 0.35); %linspace(xlim(1), xlim(2))
% lin_reg = polyfit(PDrNORMm1, PDrNORMSTN,1);
% yplot1 = polyval(lin_reg, xplot1);
% plot(xplot1,yplot1,'color','black');
% title('M1vsSTN in DYS, active');
%% Plotting M1 vs STN in Beta band
%Specifically, plotting the %Power in Beta in M1 vs the 
%Power in Beta in STN

DYSbetarest = [];
counterDYS = 0;

for i=1:numDYS
    filename = DYSdir(i).name;
    load(filename);
    if isnan(order(1,3));
        continue;
    end 
    DYSbetarest(i,1) = (subfreq(2,2,1)+subfreq(3,2,1));%#ok<AGROW> 
    DYSbetarest(i,2) = (subfreq(2,2,3)+subfreq(3,2,3));%#ok<AGROW> 
    DYSbetaactive(i,1) = (subfreq(2,4,1)+subfreq(3,4,1));%#ok<AGROW> 
    DYSbetaactive(i,2) = (subfreq(2,4,3)+subfreq(3,4,3));%#ok<AGROW>
    counterDYS = counterDYS+1;
end

DYSbetarest=100*DYSbetarest;
DYSbetaactive=100*DYSbetaactive;

[R,P] = corrcoef(DYSbetarest);
R = R(2,1);
P = P(2,1);

[Ra Pa] = corrcoef(DYSbetaactive);
Ra = Ra(2,1);
Pa = Pa(2,1);

h4=figure;
h2=subplot(2,1,1);
scatter(DYSbetarest(:,1), DYSbetarest(:,2),'MarkerEdgeColor','b')
axis([0 100 0 100]);
title('DYS Beta %Power ARM');
xlabel({'M1'});
ylabel({'STN'});
legend(['r^2 = ' num2str(R) '; p = ' num2str(P)], 'Location', 'EastOutside'); 
% hold on;
% xplot1 = linspace(0, 100); %linspace(xlim(1), xlim(2))
% lin_reg = polyfit(PDrNORMm1, PDrNORMSTN,1);
% yplot1 = polyval(lin_reg, xplot1);
% plot(xplot1,yplot1,'color','black');

h2=subplot(2,1,2);
scatter(DYSbetaactive(:,1), DYSbetaactive(:,2),'MarkerEdgeColor','r')
axis([0 100 0 100]);
legend(['r^2 = ' num2str(Ra) '; p = ' num2str(Pa)], 'Location', 'EastOutside'); 
xlabel({'M1'});
ylabel({'STN'});
% hold on;
% xplot1 = linspace(0, 100); %linspace(xlim(1), xlim(2))
% lin_reg = polyfit(PDrNORMm1, PDrNORMSTN,1);
% yplot1 = polyval(lin_reg, xplot1);
% plot(xplot1,yplot1,'color','black');