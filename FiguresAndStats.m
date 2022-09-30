%%%%%%% Figures and Stats Code 
%%%%%% Saccade Directions Head Tilt Manuscript 
%%%%%% September/October 2022

%% Initialize variables and load .mat files
% This script should live in a folder, "projectFolder", where another folder
% called "FilesForManuscript" also lives. This "FilesForManuscript" folder
% will contain many .mat files that can be used to reproduce the figures in
% our eNeuro publication. 

projectFolder = 'C:\Users\stephanie_reeves\UC Berkeley\OMlab - OM-lab-share\Projects\SaccadeDirectionsHeadTilt\Code';
SubjList = ["0101", "0122", "0123", "0124", "0125", "0126", "0127", "0128", "0129", "0166", "0167", "0169", "0172", "0173"];
DEBUG = 0; % Set this to 1 to see plots that do not appear in manuscript 

% Load a file that contains a few variables needed for figures
load(fullfile(projectFolder,'VarsForFigures'))

% Load saccade table and trial table 
load(fullfile(projectFolder,'saccade'))
load(fullfile(projectFolder,'trialTable'))

%% Head orientation tilts / numbers for manuscript 
% Get the confidence intervals for the paper
pLo = 0.05/2; %0.05 bc 95% conf intervals
pUp = 1 - 0.05/2;
crit = tinv([pLo pUp], 13); %13 bc 14-1 is the dof 

se1 = std(tab.MedianHeadTiltRightAmt,'omitnan') / sqrt(length(tab.MedianHeadTiltRightAmt));
lowerCI_tiltright = mean(tab.MedianHeadTiltRightAmt,'omitnan') + crit(1)*se1;
upperCI_tiltright = mean(tab.MedianHeadTiltRightAmt,'omitnan') + crit(2)*se1;

se2 = std(tab.MedianHeadTiltLeftAmt,'omitnan') / sqrt(length(tab.MedianHeadTiltLeftAmt));
lowerCI_tiltleft = mean(tab.MedianHeadTiltLeftAmt,'omitnan') + crit(1)*se2;
upperCI_tiltleft = mean(tab.MedianHeadTiltLeftAmt,'omitnan') + crit(2)*se2;

se3 = std(tab.MedianHeadTiltUprightAmt,'omitnan') / sqrt(length(tab.MedianHeadTiltUprightAmt));
lowerCI_upright = mean(tab.MedianHeadTiltUprightAmt,'omitnan') + crit(1)*se3;
upperCI_upright = mean(tab.MedianHeadTiltUprightAmt,'omitnan') + crit(2)*se3;


%% Main Seq figure for paper
% Load a file to look at raw traces of one subject
binsize = 10;
binedges = [0:binsize:360]/180*pi;
bincenters = [-binsize/2:binsize:360]/180*pi;
use = polarhistogram(saccade.Direction(find(saccade.TrialNumber > 0)),bincenters,'Normalization','probability');
useThis = use.Values; close;

figure('Position', [415 705 1145 633])
subplot(2,3,1)
scatter(saccade.Amplitude(find(saccade.TrialNumber > 0)), saccade.PeakSpeed(find(saccade.TrialNumber > 0)),1,'MarkerEdgeColor',[.3 .3 .3]);
xlabel('Amplitude (deg)');
ylabel('Peak Velocity (deg/s)');
ylim([0 800]);
xlim([0 60]);
set(gca,'FontSize',11)

subplot(2,3,2)
histAmp = histc(saccade.Amplitude(find(saccade.TrialNumber > 0)), [0:1:40]);
b = bar([0:1:40],histAmp,'histc');
b.FaceColor = [.6 .6 .6];
xlabel('Amplitude (deg)');
ylabel('N saccades');
ylim([0 13000]);
set(gca,'FontSize',11)

subplot(2,3,3)
%histDur = histc(saccade.DurationMs, [0:1:100]);
%bar([0:1:100],histDur,'histc');
hist = histogram(saccade.DurationMs(find(saccade.TrialNumber > 0)));
hist.BinWidth = 8;
hist.FaceColor = [.6 .6 .6]
xlim([0 500]);
ylim([0 11000]);
xlabel('Duration (ms)');
ylabel('N saccades');
set(gca,'FontSize',11)

st = 181402 
en = 182000 % 183278 is the whole trial
st_2 = 171153
en_2 = 173028
subplot(2,3,[4 5 6])
plot(samplesDataTable0169.Time(st:en)-samplesDataTable0169.Time(st), samplesDataTable0169.RightX(st:en),'Color','k') %nanmedfilt is a function that Jorge wrote for Arume and works by taking the sliding median of a certain number of samples, in this case 2000 samples
hold
plot(samplesDataTable0169.Time(st:en)-samplesDataTable0169.Time(st), samplesDataTable0169.RightY(st:en),'Color',[.6 .6 .6]) %nanmedfilt is a function that Jorge wrote for Arume and works by taking the sliding median of a certain number of samples, in this case 2000 samples
ylim([-30 30])
xlim([0 5])
xlabel('Time (s)')
ylabel('Position (deg)')
legend('Horizontal','Vertical')
set(gca,'FontSize',11)

%% Stats for manuscript
% Are the distributions for fractals in head-referenced coordinates different? 
mean(tab.xcorrKDE_HeadUpHeadRight)
mean(tab.xcorrKDE_HeadUpHeadLeft)

se1 = std(tab.xcorrKDE_HeadUpHeadRight,'omitnan') / sqrt(length(tab.xcorrKDE_HeadUpHeadRight));
lowerCI_HUHR = mean(tab.xcorrKDE_HeadUpHeadRight,'omitnan') + crit(1)*se1;
upperCI_HUHR = mean(tab.xcorrKDE_HeadUpHeadRight,'omitnan') + crit(2)*se1;

se2 = std(tab.xcorrKDE_HeadUpHeadLeft,'omitnan') / sqrt(length(tab.xcorrKDE_HeadUpHeadLeft));
lowerCI_HUHL = mean(tab.xcorrKDE_HeadUpHeadLeft,'omitnan') + crit(1)*se2;
upperCI_HUHL = mean(tab.xcorrKDE_HeadUpHeadLeft,'omitnan') + crit(2)*se2;

[h,p,ci,stats] = ttest(tab.xcorrKDE_HeadUpHeadLeft,tab.xcorrKDE_HeadUpHeadRight)



% Paired t-test to determine whether diff from zero
[h,p,ci,stats] = ttest(tab.crossCorrFractalsWorldRef_HUHL,tab.crossCorrFractalsWorldRef_HUHR)

% Get means and CIs
mean(tab.crossCorrFractalsWorldRef_HUHR)
mean(tab.crossCorrFractalsWorldRef_HUHL)
se1 = std(tab.crossCorrFractalsWorldRef_HUHR,'omitnan') / sqrt(length(tab.crossCorrFractalsWorldRef_HUHR));
lowerCI_HUHR = mean(tab.crossCorrFractalsWorldRef_HUHR,'omitnan') + crit(1)*se1;
upperCI_HUHR = mean(tab.crossCorrFractalsWorldRef_HUHR,'omitnan') + crit(2)*se1;
se2 = std(tab.crossCorrFractalsWorldRef_HUHL,'omitnan') / sqrt(length(tab.crossCorrFractalsWorldRef_HUHL));
lowerCI_HUHL = mean(tab.crossCorrFractalsWorldRef_HUHL,'omitnan') + crit(1)*se2;
upperCI_HUHL = mean(tab.crossCorrFractalsWorldRef_HUHL,'omitnan') + crit(2)*se2;

% Cohens d 
cohens = (mean(tab.crossCorrFractalsWorldRef_HUHR) - mean(tab.crossCorrFractalsWorldRef_HUHL)) / (std(tab.crossCorrFractalsWorldRef_HUHR - tab.crossCorrFractalsWorldRef_HUHL))


%% Figure (single subject) that shows circular KDE and circular cross correlation procedure 
clear vfEstimate vfEstimate2 hist1Rep hist2Rep maxlag lags x i
h1 = saccade.Direction(find(saccade.Subject == "0129" & saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.ImageTilt == -30 & saccade.TrialNumber > 0));
h2 = saccade.Direction(find(saccade.Subject == "0129" & saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.ImageTilt == 0 & saccade.TrialNumber > 0));

% Get estimates
[vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1);
[vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1);

% Do the circular cross correlation (circular bc we do the repeating)
hist1Rep = [vfEstimate vfEstimate vfEstimate];
hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
[x,lags] = xcorr(hist2Rep, hist1Rep, 450); 
[~,i] = max(x);
maxlag = lags(i)/10;

% Plot 
h = figure('Position',1.0e+03 * [0.3137    0.4637    1.2820    0.3587])
subplot(1,3,1); 
polarplot([0:.1:360]/180*pi, vfEstimate,'LineWidth',1, 'Color', [.9 .2 .2]); hold on; 
polarplot([0:.1:360]/180*pi, vfEstimate2,'LineWidth',1, 'Color', 'k');
set(gca,'fontsize',14)

subplot(1,3,2); 
plot(vfEstimate,'Color',[.9 .2 .2]); hold on; 
plot(vfEstimate2,'k');
xticks([0 900 1800 2700 3600])
xticklabels({'0','90','180','270','360'})
xlabel('Angle (deg)')
ylabel('Frequency')
set(gca,'fontsize',14)
xlim([0 3600])

subplot(1,3,3); 
plot(lags, x,'Color',[.9 .2 .2]); 
line([maxlag*10 maxlag*10],get(gca,'ylim'),'Color',[.9 .2 .2]); 
text(maxlag, max(x), sprintf('Max at %.2f deg',maxlag), 'VerticalAlignment','bottom');
xticks([-600 -400 -200 0 200 400 600])
xticklabels({'-60','-40','-20','0','20','40','60'})
xlabel('Lags (deg)')
ylabel('Correlation')
set(gca,'fontsize',14)

%% Saccade polar distributions for fractals shown both relative to the world and relative to the head
% Figure for manuscript
mean(out.median_HeadRoll(find(out.HeadTilt == "Left")));
mean(out.median_HeadRoll(find(out.HeadTilt == "Right")));

test1 = []; test2 = []; test3 = []; test4 = []; test5 = []; test6 = [];

binedges = [0:.1:360]/180*pi;

for i = 1:length(SubjList)
    h1 = saccade.DirAccountingForHeadTilt(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Left" & saccade.ImageType == "fractals" & saccade.TrialNumber > 0));
    [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test1 = [test1; vfEstimate];
    
    h2 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Upright" & saccade.ImageType == "fractals" & saccade.TrialNumber > 0));
    [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test2 = [test2; vfEstimate2];
    
    h3 = saccade.DirAccountingForHeadTilt(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Right" & saccade.ImageType == "fractals" & saccade.TrialNumber > 0)); 
    [vfEstimate3] = circ_ksdensity(h3,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test3 = [test3; vfEstimate3];
    
    h4 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Left" & saccade.ImageType == "fractals" & saccade.TrialNumber > 0));
    [vfEstimate4] = circ_ksdensity(h4,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test4 = [test4; vfEstimate4];
    
    h5 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Upright" & saccade.ImageType == "fractals" & saccade.TrialNumber > 0));
    [vfEstimate5] = circ_ksdensity(h5,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test5 = [test5; vfEstimate5];
    
    h6 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Right" & saccade.ImageType == "fractals" & saccade.TrialNumber > 0));
    [vfEstimate6] = circ_ksdensity(h6,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test6 = [test6; vfEstimate6];
    close
end

% Get the confidence intervals
pLo = 0.05/2; %0.05 bc 95% conf intervals
pUp = 1 - 0.05/2;
crit = tinv([pLo pUp], 13); %13 bc 14-1 is the dof 

se1 = std(test1,'omitnan') / sqrt(size(test1,1));
lowerCI_1 = mean(test1,'omitnan') + crit(1)*se1;
upperCI_1 = mean(test1,'omitnan') + crit(2)*se1;

se2 = std(test2,'omitnan') / sqrt(size(test2,1));
lowerCI_2 = mean(test2,'omitnan') + crit(1)*se2;
upperCI_2 = mean(test2,'omitnan') + crit(2)*se2;

se3 = std(test3,'omitnan') / sqrt(size(test3,1));
lowerCI_3 = mean(test3,'omitnan') + crit(1)*se3;
upperCI_3 = mean(test3,'omitnan') + crit(2)*se3;

se4 = std(test4,'omitnan') / sqrt(size(test4,1));
lowerCI_4 = mean(test4,'omitnan') + crit(1)*se4;
upperCI_4 = mean(test4,'omitnan') + crit(2)*se4;

se5 = std(test5,'omitnan') / sqrt(size(test5,1));
lowerCI_5 = mean(test5,'omitnan') + crit(1)*se5;
upperCI_5 = mean(test5,'omitnan') + crit(2)*se5;

se6 = std(test6,'omitnan') / sqrt(size(test6,1));
lowerCI_6 = mean(test6,'omitnan') + crit(1)*se6;
upperCI_6 = mean(test6,'omitnan') + crit(2)*se6;


h = figure('Position',1.0e+03 * [0.0010  0.0410  2.5600  1.3273])
subplot(2,3,1)
polarplot([0:.1:360]/180*pi, mean(test1),'LineWidth',2, 'Color', 'k'); hold on; ax_pol = gca;
polarplot([0; 180]*pi/180, [1 1], 'LineWidth', 2, 'Color', [.18 .78 .92]) %blue reference line
polarplot([26.7; 206.7]*pi/180, [1 1], 'LineWidth', 2, 'Color', [1 .66 .19]) %orange reference line
rlim([0 0.35])
%ax1.FontSize = 22;
%error = std(test1);
%lowerBound = plotme1-error;
%upperBound = plotme1+error;
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_1);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_1));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; 
set(ax_cart,'visible','off');
ax_pol.FontSize = 16;
nums = ax_pol.Position
ax_pol.Position = [nums(1) nums(2) nums(3)+.04 nums(4)] %shift subplot to the right
ax_cart.Position = [nums(1) nums(2) nums(3)+.04 nums(4)]

subplot(2,3,2)
polarplot([0:.1:360]/180*pi, mean(test2),'LineWidth',2, 'Color', 'k'); hold on; ax_pol = gca;
hold on
polarplot([0; 180]*pi/180, [1 1], 'LineWidth', 2, 'LineStyle', ':', 'Color', [.18 .78 .92]) %blue
polarplot([0; 180]*pi/180, [1 1], 'LineWidth', 2, 'LineStyle', '--', 'Color', [1 .66 .19]) %orange
rlim([0 0.35])
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_2);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_2));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; set(ax_cart,'visible','off')
ax_pol.FontSize = 16;

subplot(2,3,3)
polarplot([0:.1:360]/180*pi, mean(test3),'LineWidth',2, 'Color', 'k'); hold on; ax_pol = gca;
hold on
polarplot([0; 180]*pi/180, [1 1], 'LineWidth', 2, 'Color', [.18 .78 .92]) %blue
polarplot([152.9; 332.9]*pi/180, [1 1], 'LineWidth', 2, 'Color', [1 .66 .19]) %orange
rlim([0 0.35])
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_3);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_3));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; set(ax_cart,'visible','off')
ax_pol.FontSize = 16;
nums = ax_pol.Position
ax_pol.Position = [nums(1) nums(2) nums(3)-.04 nums(4)] %shift subplot to the left
ax_cart.Position = [nums(1) nums(2) nums(3)-.04 nums(4)]

subplot(2,3,4)
polarplot([0:.1:360]/180*pi, mean(test4),'LineWidth',2, 'Color', 'k'); hold on; ax_pol = gca;
hold on
polarplot([152.9; 332.9]*pi/180, [1 1], 'LineWidth', 2, 'Color', [.18 .78 .92]) %blue
polarplot([0; 180]*pi/180, [1 1], 'LineWidth', 2, 'Color', [1 .66 .19]) %orange
rlim([0 0.35])
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_4);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_4));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; set(ax_cart,'visible','off')
ax_pol.FontSize = 16;
nums = ax_pol.Position;
ax_pol.Position = [nums(1) nums(2) nums(3)+.04 nums(4)]; %shift subplot to the right
ax_cart.Position = [nums(1) nums(2) nums(3)+.04 nums(4)];

subplot(2,3,5)
polarplot([0:.1:360]/180*pi, mean(test5),'LineWidth',2, 'Color', 'k'); hold on; ax_pol = gca;
hold on
polarplot([0; 180]*pi/180, [1 1], 'LineWidth', 2, 'LineStyle', ':', 'Color', [.18 .78 .92]) %blue
polarplot([0; 180]*pi/180, [1 1], 'LineWidth', 2, 'LineStyle', '--', 'Color', [1 .66 .19]) %orange
rlim([0 0.35])
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_5);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_5));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; set(ax_cart,'visible','off')
ax_pol.FontSize = 16;

subplot(2,3,6)
polarplot([0:.1:360]/180*pi, mean(test6),'LineWidth',2, 'Color', 'k'); hold on; ax_pol = gca;
hold on
polarplot([26.7; 206.7]*pi/180, [1 1], 'LineWidth', 2, 'Color', [.18 .78 .92]) %blue
polarplot([0; 180]*pi/180, [1 1], 'LineWidth', 2, 'Color', [1 .66 .19]) %orange
rlim([0 0.35])
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_6);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_6));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; set(ax_cart,'visible','off')
ax_pol.FontSize = 16;
nums = ax_pol.Position;
ax_pol.Position = [nums(1) nums(2) nums(3)-.04 nums(4)]; %shift subplot to the left
ax_cart.Position = [nums(1) nums(2) nums(3)-.04 nums(4)];

%%%%%%% When I submitted this manuscript originally, I just saved figures as SVG 
%%%%%%% and was able to export with transparency as vector format. 

%%%%%% For the second submission (eNeuro), you HAVE to include this line
%%%%%% before saving as normal as SVF.
set(gcf,'renderer','Painters') 
%%%%%% ^^ that is the line to include!!! 


%% MANUSCRIPT FIGURE -- Quantifying the Offset for Fractals
% Need to iterate through for accurate figure w/ shading
binedges = [0:.1:360]/180*pi;
test1 = [];
test2 = [];
test3 = [];

% Get the vals for each subject and stitch them together
for i = 1:length(SubjList)
    h1 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Left" & saccade.ImageType == "fractals" & saccade.TrialNumber > 0));
    [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test1 = [test1; vfEstimate];
    
    h2 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Right" & saccade.ImageType == "fractals" & saccade.TrialNumber > 0));
    [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test2 = [test2; vfEstimate2];
    
    h3 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Upright" & saccade.ImageType == "fractals" & saccade.TrialNumber > 0));
    [vfEstimate3] = circ_ksdensity(h3,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test3 = [test3; vfEstimate3];
    close
end

% Get the confidence intervals
pLo = 0.05/2; %0.05 bc 95% conf intervals
pUp = 1 - 0.05/2;
crit = tinv([pLo pUp], 13); %13 bc 14-1 is the dof 

se1 = std(test1,'omitnan') / sqrt(size(test1,1));
lowerCI_1 = mean(test1,'omitnan') + crit(1)*se1;
upperCI_1 = mean(test1,'omitnan') + crit(2)*se1;

se2 = std(test2,'omitnan') / sqrt(size(test2,1));
lowerCI_2 = mean(test2,'omitnan') + crit(1)*se2;
upperCI_2 = mean(test2,'omitnan') + crit(2)*se2;

se3 = std(test3,'omitnan') / sqrt(size(test3,1));
lowerCI_3 = mean(test3,'omitnan') + crit(1)*se3;
upperCI_3 = mean(test3,'omitnan') + crit(2)*se3;

% Make figure for manuscript
figure('Position',1.0e+03 * [0.0100  0.0610  1.1243  0.9867])
sp1 = subplot(2,10,1:5)
sp1.Position = sp1.Position - [0 0.06 0 0];
polarplot(binedges,mean(test3),'LineWidth',2, 'Color', 'k') ; hold on; ax_pol = gca;
polarplot(binedges,mean(test2),'LineWidth',2, 'Color', [.9 .2 .2]) %red
rlim([0 0.35])
leg = legend('Head Upright','Head Right')
leg.Location = 'south';
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_2);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_2));
fill([xl,xh],[yl,yh],[.9 .2 .2],'FaceAlpha',0.5,'EdgeAlpha',0)
hold
[xl,yl] = pol2cart(binedges,lowerCI_3);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_3));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; set(ax_cart,'visible','off')

sp2 = subplot(2,10,6:10)
sp2.Position = sp2.Position - [0 0.06 0 0];
polarplot(binedges,mean(test3),'LineWidth',2, 'Color', 'k'); hold on; ax_pol = gca; 
polarplot(binedges,mean(test1),'LineWidth',2, 'Color', [.1 .5 .8]) %blue
rlim([0 0.35])
leg = legend('Head Upright','Head Left')
leg.Location = 'south';
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_1);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_1));
fill([xl,xh],[yl,yh],[.1 .5 .8],'FaceAlpha',0.5,'EdgeAlpha',0)
hold
[xl,yl] = pol2cart(binedges,lowerCI_3);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_3));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; set(ax_cart,'visible','off')

sp3 = subplot(2,10,11:15)
bar([1.2:length(SubjList)+.2],bootstrapped_right.Avg,'FaceColor',[.9 .2 .2],'EdgeColor',[.9 .2 .2], 'LineWidth',.5); hold on;
bar([1:length(SubjList)],bootstrapped_left.Avg,'FaceColor',[.1 .5 .8],'EdgeColor',[.1 .5 .8], 'LineWidth',.5); hold on;
for subj = 1:length(SubjList)
    line([subj+.2 subj+.2],[bootstrapped_right.CI_low(subj) bootstrapped_right.CI_high(subj) ],'Color','k','LineWidth', .75)
    line([subj subj],[bootstrapped_left.CI_low(subj) bootstrapped_left.CI_high(subj) ],'Color','k','LineWidth', .75)
end
xticks(1:length(SubjList))
xticklabels(bootstrapped_right.SubjNum'); %set(gca,'XTickLabelRotation',80);
xlabel('Individual Subjects')
ylabel('Direction Distribution Displacement (deg)')
ylim([-34 34])
set(gca,'FontSize',12);


sp4 = subplot(2,10,16:17)
sp4.Position = sp4.Position - [0.02 0 0 0]
p1x = 1;
p2x = 2;
p1y = mean(tab.xcorrKDE_HeadUpHeadRight);
p2y = mean(tab.xcorrKDE_HeadUpHeadLeft);
% Get the confidence intervals
se1 = std(tab.xcorrKDE_HeadUpHeadRight,'omitnan') / sqrt(length(tab.xcorrKDE_HeadUpHeadRight));
ci1 = mean(tab.xcorrKDE_HeadUpHeadRight,'omitnan') + crit*se1;
se2 = std(tab.xcorrKDE_HeadUpHeadLeft,'omitnan') / sqrt(length(tab.xcorrKDE_HeadUpHeadLeft));
ci2 = mean(tab.xcorrKDE_HeadUpHeadLeft,'omitnan') + crit*se2;
p1 = bar(p1x, p1y,'FaceColor',[1 1 1],'EdgeColor',[.9 .2 .2], 'LineWidth',2) % edge color used to be .3 .3 .3
hold on
p2 = bar(p2x, p2y,'FaceColor',[1 1 1],'EdgeColor',[.1 .5 .8], 'LineWidth',2) % edge color used to be .3 .3 .3
for i = 1:length(tab.xcorrKDE_HeadUpHeadRight)
    scatter(p1x, tab.xcorrKDE_HeadUpHeadRight(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
for i = 1:length(tab.xcorrKDE_HeadUpHeadLeft)
    scatter(p2x, tab.xcorrKDE_HeadUpHeadLeft(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
line([1 1],[ci1(1) ci1(2)],'Color','k','LineWidth', .75)
line([2 2],[ci2(1) ci2(2)],'Color','k','LineWidth', .75)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',12,'XTickLabelRotation',0,'yticklabels',[])
xticks(1:2)
xticklabels({'Right','Left'})
xlabel('Head Tilts')
ylim([-34 34])

sp5 = subplot(2,10,18:20)
RFindex_frac = (mean([tab.xcorrKDE_HeadUpHeadLeft -tab.xcorrKDE_HeadUpHeadRight],2))./(mean([-MedianHeadTiltLeftAmt MedianHeadTiltRightAmt],2));
bp = boxplot(RFindex_frac, 'Colors','k','Symbol','+k') %the +k means use the symbol + for the outlier and use the color k for the outlier
set(bp,'LineWidth',1)
hold
for i = 1:length(RFindex_frac)
    scatter(1, RFindex_frac(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
set(gca,'FontSize',12);
ylim([-0.14 1.1])
xticklabels({' '})
ylabel('Reference Frame Index')
txt = {'     Head','\leftarrow Orientation', '     Hypothesis'};
text(.5, 0, txt, 'FontSize',11)
txt2 = {'\leftarrow World Orientation','     Hypothesis'};
text(.5, 0.975, txt2, 'FontSize',11)

%%%%%% For the second submission (eNeuro), you HAVE to include this line
%%%%%% before saving as normal as SVF.
set(gcf,'renderer','Painters') 

mean(RFindex_frac)
[h,p,ci,stats] = ttest(RFindex_frac)

% Cohens d 
cohens = (mean(RFindex_frac) - 0) / (std(RFindex_frac))
%cohens = (mean(xcorrKDE_HeadUpHeadLeft) - mean(xcorrKDE_HeadUpHeadRight)) / (std(xcorrKDE_HeadUpHeadLeft) - std(xcorrKDE_HeadUpHeadRight))


%% MANUSCRIPT FIGURE -- Torsion for JUST fractals, including torsion session traces for one subject
% Load this file for one subject so that you can show the raw torsion traces
Left_start = find(samplesDataTable0173.TrialNumber == 1,1); % find the row where the first instance of trial 1 occurred
Left_end = find(samplesDataTable0173.TrialNumber == 60,1,'last');
Upright_start = find(samplesDataTable0173.TrialNumber == 61,1);
Upright_end = find(samplesDataTable0173.TrialNumber == 120,1,'last');
Right_start = find(samplesDataTable0173.TrialNumber == 121,1);
Right_end = find(samplesDataTable0173.TrialNumber == 180,1,'last');

% Calculate other info
fractalsHeadLeftTorsion = statarray4.mean_median_Torsion(find(statarray4.HeadTilt == "Left" & statarray4.ImageType == "fractals" & statarray4.ImageTilt == 0)); %the image tilt here isnt doing anything.... just a way to get from 3 repeating vals per subject to 1 single val
fractalsHeadRightTorsion = statarray4.mean_median_Torsion(find(statarray4.HeadTilt == "Right" & statarray4.ImageType == "fractals" & statarray4.ImageTilt == 0));
fractalsHeadUprightTorsion = statarray4.mean_median_Torsion(find(statarray4.HeadTilt == "Upright" & statarray4.ImageType == "fractals" & statarray4.ImageTilt == 0));
fractalsUMinusR_Torsion = fractalsHeadRightTorsion - fractalsHeadUprightTorsion;
fractalsUMinusL_Torsion = fractalsHeadLeftTorsion - fractalsHeadUprightTorsion;
p1x = 1;
p1y = mean(fractalsUMinusR_Torsion,'omitnan')
p2x = 2;
p2y = mean(fractalsUMinusL_Torsion,'omitnan')

% Get the confidence intervals
pLo = 0.05/2; %0.05 bc 95% conf intervals
pUp = 1 - 0.05/2;
crit = tinv([pLo pUp], 13); %13 bc 14-1 is the dof 
se1 = std(fractalsUMinusR_Torsion,'omitnan') / sqrt(length(fractalsUMinusR_Torsion));
ci1 = mean(fractalsUMinusR_Torsion,'omitnan') + crit*se1;
se2 = std(fractalsUMinusL_Torsion,'omitnan') / sqrt(length(fractalsUMinusL_Torsion));
ci2 = mean(fractalsUMinusL_Torsion,'omitnan') + crit*se2;
% stderr2 = std(fractalsUMinusL, 'omitnan') / sqrt(length(fractalsUMinusL));

figure('Position', 1.0e+03 * [0.0797  0.1283  1.2567  0.6600])
sp123 = subplot(2,3,[1 2 3])
plot(samplesDataTable0173.Time(Left_start:Left_end)-samplesDataTable0173.Time(Left_start), nanmedfilt( samplesDataTable0173.RightT(Left_start:Left_end),3000),'Color',[.1 .5 .8]) %nanmedfilt is a function that Jorge wrote for Arume and works by taking the sliding median of a certain number of samples, in this case 2000 samples
hold on
plot(samplesDataTable0173.Time(Right_start:Right_end)-samplesDataTable0173.Time(Right_start), nanmedfilt( samplesDataTable0173.RightT(Right_start:Right_end),2000),'Color',[.9 .2 .2])
plot(samplesDataTable0173.Time(Upright_start:Upright_end)-samplesDataTable0173.Time(Upright_start), nanmedfilt( samplesDataTable0173.RightT(Upright_start:Upright_end),2000),'Color','k')
ylim([-10 10])
xlim([0 1000])
ylabel('Ocular Counter Roll (deg)','FontSize',12)
xlabel('Time (s)','FontSize',12)
legend('Head Left','Head Right','Head Upright')

sp4 = subplot(2,3,4)
p1 = bar(p1x, p1y, 'FaceColor',[1 1 1],'EdgeColor',[.9 .2 .2],'LineWidth', 2)
hold on
p2 = bar(p2x, p2y, 'FaceColor',[1 1 1],'EdgeColor',[.1 .5 .8],'LineWidth', 2) 
for i = 1:length(fractalsUMinusR_Torsion)
    scatter(p1x, fractalsUMinusR_Torsion(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
for i = 1:length(fractalsUMinusL_Torsion)
    scatter(p2x, fractalsUMinusL_Torsion(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
line([1 1],[ci1(1) ci1(2)],'Color','k','LineWidth', .7)
line([2 2],[ci2(1) ci2(2)],'Color','k','LineWidth', .7)
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'fontsize',12)
xticks(1:2)
xticklabels({'Right','Left'})
xlabel('Head Tilts','FontSize',12);
ylabel('Change in Median OCR (deg)','FontSize',12)
ylim([-11 11])

sp5 = subplot(2,3,5)
scatter(-tab.xcorrKDE_HeadUpHeadRight, -fractalsUMinusR_Torsion, 'MarkerEdgeColor',[.9 .2 .2], 'LineWidth', 1)
hold
scatter(tab.xcorrKDE_HeadUpHeadLeft, fractalsUMinusL_Torsion, 'MarkerEdgeColor',[.1 .5 .8], 'LineWidth', 1)
xlabel('Direction Distribution Displacement (deg)','FontSize',12);
ylabel('Change in Median OCR (deg)','FontSize',12);
xlim([-5 17]);
ylim([-2 11]);
box on;

sp6 = subplot(2,3,6)
sp6.Position = sp6.Position + [-0.0227 0 0.0227 0];
RFI = (mean([-tab.xcorrKDE_HeadUpHeadRight tab.xcorrKDE_HeadUpHeadLeft],2))./(mean([fractalsUMinusL_Torsion -fractalsUMinusR_Torsion],2));
bp = boxplot(RFI, 'Colors','k','Symbol','+k') %the +k means use the symbol + for the outlier and use the color k for the outlier
hold
scatter(1,RFI,20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
set(bp,'LineWidth',1.5)
xticklabels({' '})
ylabel('Reference Frame Index','FontSize',12)
txt = {'\leftarrow Head'};
text(.5, -0.025, txt, 'FontSize',11)
txt2 = {'\leftarrow Retinal'};
text(.5, 0.975, txt2, 'FontSize',11)
ylim([-.25 1.6])

mean(fractalsUMinusR_Torsion,'omitnan')
mean(fractalsUMinusL_Torsion,'omitnan')

se1 = std(fractalsUMinusR_Torsion,'omitnan') / sqrt(length(fractalsUMinusR_Torsion));
lowerCI_torsionR = mean(fractalsUMinusR_Torsion,'omitnan') + crit(1)*se1;
upperCI_torsionR = mean(fractalsUMinusR_Torsion,'omitnan') + crit(2)*se1;

se2 = std(fractalsUMinusL_Torsion,'omitnan') / sqrt(length(fractalsUMinusL_Torsion));
lowerCI_torsionL = mean(fractalsUMinusL_Torsion,'omitnan') + crit(1)*se2;
upperCI_torsionL = mean(fractalsUMinusL_Torsion,'omitnan') + crit(2)*se2;

mean(RFI,'omitnan')
[h,p,ci,stats] = ttest(RFI)
[h,p,ci,stats] = ttest(1-RFI)

% Cohen's D
cohens = mean(RFI,'omitnan') / std(RFI,'omitnan')
cohens = mean(1-RFI,'omitnan') / std(1-RFI,'omitnan')

% Do correlation of scatter, taking out nans first 
idx=find(isnan(fractalsUMinusR_Torsion));
fractalsUMinusR_Torsion(idx) = [];
xcorrKDE_HeadUpHeadRight(idx) = [];

idx=find(isnan(fractalsUMinusL_Torsion));
fractalsUMinusL_Torsion(idx) = [];
xcorrKDE_HeadUpHeadLeft(idx) = [];

newT = [fractalsUMinusL_Torsion; -fractalsUMinusR_Torsion];
newxCorr=[xcorrKDE_HeadUpHeadLeft';-xcorrKDE_HeadUpHeadRight'];
[rho,pval] = corr(newT,newxCorr)



%% Earth upright scenes at diff head tilt -- figure for manuscript
% Need to iterate through the shading
binedges = [0:.1:360]/180*pi;
test1 = []; test2 = []; test3 = [];

% Get the vals for each subject and stitch them together
for i = 1:length(SubjList)
    h1 = saccade.DirAccountingForHeadTilt(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Left" & saccade.ImageType == "scenes" & saccade.ImageTilt == 30 & saccade.TrialNumber > 0));
    [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test1 = [test1; vfEstimate];
    
    h2 = saccade.DirAccountingForHeadTilt(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Right" & saccade.ImageType == "scenes" & saccade.ImageTilt == -30 & saccade.TrialNumber > 0));
    [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test2 = [test2; vfEstimate2];
    
    h3 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.ImageTilt == 0 & saccade.TrialNumber > 0));
    [vfEstimate3] = circ_ksdensity(h3,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test3 = [test3; vfEstimate3];
    close
end

% Get the confidence intervals
se1 = std(test1,'omitnan') / sqrt(size(test1,1));
lowerCI_1 = mean(test1,'omitnan') + crit(1)*se1;
upperCI_1 = mean(test1,'omitnan') + crit(2)*se1;

se2 = std(test2,'omitnan') / sqrt(size(test2,1));
lowerCI_2 = mean(test2,'omitnan') + crit(1)*se2;
upperCI_2 = mean(test2,'omitnan') + crit(2)*se2;

se3 = std(test3,'omitnan') / sqrt(size(test3,1));
lowerCI_3 = mean(test3,'omitnan') + crit(1)*se3;
upperCI_3 = mean(test3,'omitnan') + crit(2)*se3;


% Make the figure!
figure('Position',1.0e+03 * [ 0.0100    0.0610    1.0497    0.9867])

sp1 = subplot(2,10,1:5)
sp1.Position = sp1.Position - [0 0.06 0 0];
polarplot(binedges,mean(test3),'LineWidth',2, 'Color', 'k'); hold on; ax_pol = gca;
polarplot(binedges,mean(test2),'LineWidth',2, 'Color', [.9 .2 .2]) %red
rlim([0 0.4])
leg = legend('Head Upright','Head Right')
leg.Location = 'south';
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_3);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_3));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
hold
[xl,yl] = pol2cart(binedges,lowerCI_2);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_2));
fill([xl,xh],[yl,yh],[.9 .2 .2],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; set(ax_cart,'visible','off')

sp2 = subplot(2,10,6:10)
sp2.Position = sp2.Position - [0 0.06 0 0];
polarplot(binedges,mean(test3),'LineWidth',2, 'Color', 'k'); hold on; ax_pol = gca;
polarplot(binedges,mean(test1),'LineWidth',2, 'Color', [.1 .5 .8]) %blue
rlim([0 0.4])
leg = legend('Head Upright','Head Left')
leg.Location = 'south';
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_3);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_3));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
hold
[xl,yl] = pol2cart(binedges,lowerCI_1);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_1));
fill([xl,xh],[yl,yh],[.1 .5 .8],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; set(ax_cart,'visible','off')


sp3 = subplot(2,10,11:15)
bar([1:length(SubjList)],bootstrapped_30.Avg,'FaceColor',[.9 .2 .2],'EdgeColor',[.9 .2 .2], 'LineWidth',.5); hold on;
bar([1:length(SubjList)],bootstrapped_neg30.Avg,'FaceColor',[.1 .5 .8],'EdgeColor',[.1 .5 .8], 'LineWidth',.5); hold on;
for subj = 1:length(SubjList)
    line([subj subj],[bootstrapped_30.CI_low(subj) bootstrapped_30.CI_high(subj) ],'Color','k','LineWidth', .75)
    line([subj+.1 subj+.1],[bootstrapped_neg30.CI_low(subj) bootstrapped_neg30.CI_high(subj) ],'Color','k','LineWidth', .75)
end
xticks(1:length(SubjList))
xticklabels(bootstrapped_30.SubjNum'); 
yticks(-50:10:50)
yticklabels({-50:10:50})
xlabel('Individual Subjects')
ylabel('Direction Distribution Displacement (deg)')
ylim([-50 50])
set(gca,'FontSize',12);


sp4 = subplot(2,10,16:17)
sp4.Position = sp4.Position - [0.02 0 0 0]
p1x = 1;
p2x = 2;
p1y = mean(tab.crossCorrHeadUprightandRightWithImageEarthUpright);
p2y = mean(tab.crossCorrHeadUprightandLeftWithImageEarthUpright);
% Get the confidence intervals
pLo = 0.05/2; %0.05 bc 95% conf intervals
pUp = 1 - 0.05/2;
crit = tinv([pLo pUp], 13); %13 bc 14-1 is the dof 
se1 = std(tab.crossCorrHeadUprightandRightWithImageEarthUpright,'omitnan') / sqrt(length(tab.crossCorrHeadUprightandRightWithImageEarthUpright));
ci1 = mean(tab.crossCorrHeadUprightandRightWithImageEarthUpright,'omitnan') + crit*se1;
se2 = std(tab.crossCorrHeadUprightandLeftWithImageEarthUpright,'omitnan') / sqrt(length(tab.crossCorrHeadUprightandLeftWithImageEarthUpright));
ci2 = mean(tab.crossCorrHeadUprightandLeftWithImageEarthUpright,'omitnan') + crit*se2;
p1 = bar(p1x, p1y,'FaceColor',[1 1 1],'EdgeColor',[.9 .2 .2], 'LineWidth',2) 
hold on
p2 = bar(p2x, p2y,'FaceColor',[1 1 1],'EdgeColor',[.1 .5 .8], 'LineWidth',2) 
for i = 1:length(tab.crossCorrHeadUprightandRightWithImageEarthUpright)
    scatter(p1x, tab.crossCorrHeadUprightandRightWithImageEarthUpright(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
for i = 1:length(tab.crossCorrHeadUprightandLeftWithImageEarthUpright)
    scatter(p2x, tab.crossCorrHeadUprightandLeftWithImageEarthUpright(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
line([1 1],[ci1(1) ci1(2)],'Color','k','LineWidth', .75)
line([2 2],[ci2(1) ci2(2)],'Color','k','LineWidth', .75)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',12,'XTickLabelRotation',0)
xticks(1:2)
xticklabels({'Right','Left'})
set(gca,'yticklabels',[])
xlabel('Head Tilts')
ylim([-50 50])


sp5 = subplot(2,10,18:20)
RFindex_earth = (mean([-tab.crossCorrHeadUprightandLeftWithImageEarthUpright tab.crossCorrHeadUprightandRightWithImageEarthUpright],2))./mean([-MedianHeadTiltLeftAmt MedianHeadTiltRightAmt],2);
bp = boxplot(RFindex_earth, 'Colors','k','Symbol','+k') %the +k means use the symbol + for the outlier and use the color k for the outlier
set(bp,'LineWidth',1.5);
hold
for i = 1:length(RFindex_earth)
    scatter(1, RFindex_earth(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
set(gca,'FontSize',12);
xticklabels({' '});
ylabel('Reference Frame Index');
txt = {'\leftarrow Image Orientation', '     Hypothesis'};
text(.5, -0.025, txt, 'FontSize',11);
txt2 = {'\leftarrow Head Orientation','     Hypothesis'};
text(.5, 0.975, txt2, 'FontSize',11);
ylim([-0.12 1.1])

mean(RFindex_earth)
[h,p,ci,stats] = ttest(RFindex_earth)
[h,p,ci,stats] = ttest(1-RFindex_earth)

cohens = (mean(RFindex_earth)) / (std(RFindex_earth))
cohens = (mean(1-RFindex_earth)) / (std(1-RFindex_earth))


%%%%%% For the second submission (eNeuro), you HAVE to include this line
%%%%%% before saving as normal as SVF.
set(gcf,'renderer','Painters') 

% Get means and confidence intervals -- repeat from above i believe
mean(tab.crossCorrHeadUprightandRightWithImageEarthUpright)
mean(tab.crossCorrHeadUprightandLeftWithImageEarthUpright)
se1 = std(tab.crossCorrHeadUprightandRightWithImageEarthUpright,'omitnan') / sqrt(length(tab.crossCorrHeadUprightandRightWithImageEarthUpright));
lowerCI_right = mean(tab.crossCorrHeadUprightandRightWithImageEarthUpright,'omitnan') + crit(1)*se1;
upperCI_right = mean(tab.crossCorrHeadUprightandRightWithImageEarthUpright,'omitnan') + crit(2)*se1;
se2 = std(tab.crossCorrHeadUprightandLeftWithImageEarthUpright,'omitnan') / sqrt(length(tab.crossCorrHeadUprightandLeftWithImageEarthUpright));
lowerCI_left = mean(tab.crossCorrHeadUprightandLeftWithImageEarthUpright,'omitnan') + crit(1)*se2;
upperCI_left = mean(tab.crossCorrHeadUprightandLeftWithImageEarthUpright,'omitnan') + crit(2)*se2;


%% Image Tilt effect during head upright -- figure for manuscript
% Need to iterate through code in order to get accurate shading 
binedges = [0:.1:360]/180*pi;
test1 = []; test2 = []; test3 = [];

% Get the vals for each subject and stitch them together
for i = 1:length(SubjList)
    h1 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.ImageTilt == 0 & saccade.TrialNumber > 0));
    [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test1 = [test1; vfEstimate];
    
    h2 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.ImageTilt == 30 & saccade.TrialNumber > 0));
    [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test2 = [test2; vfEstimate2];
    
    h3 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.ImageTilt == -30 & saccade.TrialNumber > 0));
    [vfEstimate3] = circ_ksdensity(h3,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test3 = [test3; vfEstimate3];
    close
end

% Get the confidence intervals
se1 = std(test1,'omitnan') / sqrt(size(test1,1));
lowerCI_1 = mean(test1,'omitnan') + crit(1)*se1;
upperCI_1 = mean(test1,'omitnan') + crit(2)*se1;

se2 = std(test2,'omitnan') / sqrt(size(test2,1));
lowerCI_2 = mean(test2,'omitnan') + crit(1)*se2;
upperCI_2 = mean(test2,'omitnan') + crit(2)*se2;

se3 = std(test3,'omitnan') / sqrt(size(test3,1));
lowerCI_3 = mean(test3,'omitnan') + crit(1)*se3;
upperCI_3 = mean(test3,'omitnan') + crit(2)*se3;



figure('Position',1.0e+03 * [ 0.0100    0.0610    1.0497    0.9867])
    
sp1 = subplot(2,10,1:5)
sp1.Position = sp1.Position - [0 0.06 0 0];
polarplot(binedges,mean(test1),'LineWidth',2, 'Color', 'k') ; hold on; ax_pol = gca;
polarplot(binedges,mean(test2),'LineWidth',2, 'Color', [.5 .8 .3]) %green
rlim([0 0.4])
leg = legend('0 deg Scene Tilt','30 deg Scene Tilt')
leg.Location = 'south';
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_1);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_1));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
hold
[xl,yl] = pol2cart(binedges,lowerCI_2);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_2));
fill([xl,xh],[yl,yh],[.5 .8 .3],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; set(ax_cart,'visible','off')


sp2 = subplot(2,10,6:10)
sp2.Position = sp2.Position - [0 0.06 0 0];
polarplot(binedges,mean(test1),'LineWidth',2, 'Color', 'k'); hold on; ax_pol = gca;
polarplot(binedges,mean(test3),'LineWidth',2, 'Color', [.9 .7 .1]) %gold
rlim([0 0.4])
leg = legend('0 deg Scene Tilt','-30 deg Scene Tilt')
leg.Location = 'south';
ax_cart = axes();
ax_cart.Position = ax_pol.Position;
[xl,yl] = pol2cart(binedges,lowerCI_1);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_1));
fill([xl,xh],[yl,yh],[.6 .6 .6],'FaceAlpha',0.5,'EdgeAlpha',0)
hold
[xl,yl] = pol2cart(binedges,lowerCI_3);
[xh,yh] = pol2cart(fliplr(binedges),fliplr(upperCI_3));
fill([xl,xh],[yl,yh],[.9 .7 .1],'FaceAlpha',0.5,'EdgeAlpha',0)
xlim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
ylim(ax_cart,[-max(get(ax_pol,'RLim')),max(get(ax_pol,'RLim'))]);
axis square; set(ax_cart,'visible','off')


sp3 = subplot(2,10,11:15)
bar([1:length(SubjList)],bootstrapped_30.Avg,'FaceColor',[.5 .8 .3],'EdgeColor',[.5 .8 .3], 'LineWidth',.5); hold on;
bar([1:length(SubjList)],bootstrapped_neg30.Avg,'FaceColor',[.9 .7 .1],'EdgeColor',[.9 .7 .1], 'LineWidth',.5); hold on;
for subj = 1:length(SubjList)
    line([subj subj],[bootstrapped_30.CI_low(subj) bootstrapped_30.CI_high(subj) ],'Color','k','LineWidth', .75)
    line([subj+.1 subj+.1],[bootstrapped_neg30.CI_low(subj) bootstrapped_neg30.CI_high(subj) ],'Color','k','LineWidth', .75)
end
xticks(1:length(SubjList))
xticklabels(bootstrapped_30.SubjNum'); 
yticks(-50:10:50)
yticklabels({-50:10:50})
xlabel('Individual Subjects')
ylabel('Direction Distribution Displacement (deg)')
ylim([-50 50])
set(gca,'FontSize',12);


sp4 = subplot(2,10,16:17)
sp4.Position = sp4.Position - [0.02 0 0 0]
p1x = 1;
p2x = 2;
p1y = mean(tab.crossCorrHeadUprightImageTilt30and0);
p2y = mean(tab.crossCorrHeadUprightImageTiltneg30and0);
% Get the confidence intervals
pLo = 0.05/2; %0.05 bc 95% conf intervals
pUp = 1 - 0.05/2;
crit = tinv([pLo pUp], 13); %13 bc 14-1 is the dof 
se1 = std(tab.crossCorrHeadUprightImageTilt30and0,'omitnan') / sqrt(length(tab.crossCorrHeadUprightImageTilt30and0));
ci1 = mean(tab.crossCorrHeadUprightImageTilt30and0,'omitnan') + crit*se1;
se2 = std(tab.crossCorrHeadUprightImageTiltneg30and0,'omitnan') / sqrt(length(tab.crossCorrHeadUprightImageTiltneg30and0));
ci2 = mean(tab.crossCorrHeadUprightImageTiltneg30and0,'omitnan') + crit*se2;
p1 = bar(p1x, p1y,'FaceColor',[1 1 1],'EdgeColor',[.5 .8 .3], 'LineWidth',2) % edge color used to be .3 .3 .3
hold on
p2 = bar(p2x, p2y,'FaceColor',[1 1 1],'EdgeColor',[.9 .7 .1], 'LineWidth',2) % edge color used to be .3 .3 .3
for i = 1:length(tab.crossCorrHeadUprightImageTilt30and0)
    scatter(p1x, tab.crossCorrHeadUprightImageTilt30and0(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
for i = 1:length(tab.crossCorrHeadUprightImageTiltneg30and0)
    scatter(p2x, tab.crossCorrHeadUprightImageTiltneg30and0(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
line([1 1],[ci1(1) ci1(2)],'Color','k','LineWidth', .75)
line([2 2],[ci2(1) ci2(2)],'Color','k','LineWidth', .75)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',12,'XTickLabelRotation',0,'yticklabels',[])
xticks(1:2)
xticklabels({'30','-30'})
xlabel('Image Tilts')
ylim([-50 50])

sp5 = subplot(2,10,18:20)
RFindex_im = (mean([tab.crossCorrHeadUprightImageTilt30and0 -tab.crossCorrHeadUprightImageTiltneg30and0],2))./30;
bp = boxplot(RFindex_im, 'Colors','k','Symbol','+k') %the +k means use the symbol + for the outlier and use the color k for the outlier
set(bp,'LineWidth',1.5)
hold
for i = 1:length(RFindex_im)
    scatter(1, RFindex_im(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
set(gca,'FontSize',12);
ylim([-0.12 1.1])
xticklabels({' '})
ylabel('Reference Frame Index')
txt = {'\leftarrow Head Orientation', '     Hypothesis'};
text(.5, -0.025, txt, 'FontSize',11)
txt2 = {'\leftarrow Image Orientation','     Hypothesis'};
text(.5, 0.975, txt2, 'FontSize',11)



%%%%%% For the second submission (eNeuro), you HAVE to include this line
%%%%%% before saving as normal as SVF.
set(gcf,'renderer','Painters') 



mean(RFindex_im)
[h,p,ci,stats] = ttest(RFindex_im)
[h,p,ci,stats] = ttest(1-RFindex_im)

cohens = (mean(RFindex_im)) / (std(RFindex_im))
cohens = (mean(1-RFindex_im)) / (std(1-RFindex_im))

mean(tab.crossCorrHeadUprightImageTilt30and0)
mean(tab.crossCorrHeadUprightImageTiltneg30and0)

se1 = std(tab.crossCorrHeadUprightImageTilt30and0,'omitnan') / sqrt(length(tab.crossCorrHeadUprightImageTilt30and0));
lowerCI_30 = mean(tab.crossCorrHeadUprightImageTilt30and0,'omitnan') + crit(1)*se1;
upperCI_30 = mean(tab.crossCorrHeadUprightImageTilt30and0,'omitnan') + crit(2)*se1;

se2 = std(tab.crossCorrHeadUprightImageTiltneg30and0,'omitnan') / sqrt(length(tab.crossCorrHeadUprightImageTiltneg30and0));
lowerCI_neg30 = mean(tab.crossCorrHeadUprightImageTiltneg30and0,'omitnan') + crit(1)*se2;
upperCI_neg30 = mean(tab.crossCorrHeadUprightImageTiltneg30and0,'omitnan') + crit(2)*se2;

%% Image Tilt effect across other head tilts? Stats

% Get the RF indicies 
RFindex_upright = (mean([tab.crossCorrHeadUprightImageTilt30and0 -tab.crossCorrHeadUprightImageTiltneg30and0],2))./30;
RFindex_right = (mean([tab.crossCorrHeadRightImageTilt30and0 -tab.crossCorrHeadRightImageTiltneg30and0],2))./mean(MedianHeadTiltRightAmt);
RFindex_left = (mean([tab.crossCorrHeadLeftImageTilt30and0 -tab.crossCorrHeadLeftImageTiltneg30and0],2))./mean(-MedianHeadTiltLeftAmt);


% Statistical differences? 
forAnova = [RFindex_left RFindex_upright RFindex_right];
[p,tbl,stats] = anova1(forAnova)
results = multcompare(stats);
tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])







