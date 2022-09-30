%%%%%% Analysis Code 
%%%%%% Saccade Directions Head Tilt Manuscript 
%%%%%% September/October 2022
%%%%%% Stephanie Reeves, Otero-Millan Laboratory

% Initialize variables 
projectFolder = 'C:\Users\stephanie_reeves\UC Berkeley\OMlab - OM-lab-share\Projects\SaccadeDirectionsHeadTilt\Code';
saveHere = 'C:\Users\stephanie_reeves\UC Berkeley\OMlab - OM-lab-share\Projects\SaccadeDirectionsHeadTilt\Code\FilesForManuscript';
SubjList = ["0101", "0122", "0123", "0124", "0125", "0126", "0127", "0128", "0129", "0166", "0167", "0169", "0172", "0173"];
subjects = {'0101', '0122', '0123', '0124', '0125', '0126', '0127', '0128', '0129', '0166',  '0167',  '0169',  '0172',  '0173'};
numSubjList = [101, 122, 123, 124, 125, 126, 127, 128, 129, 166, 167, 169, 172, 173]'; % should match above line
DEBUG = 0; % Set this to 1 to see plots that do not appear in manuscript 
saveFiles = 0; % Set this to 1 if you want to save some .mat files/tables

%% Load saccade and trial table
load(fullfile(projectFolder,'saccade'))
load(fullfile(projectFolder,'trialTable'))

% Get torsion stuff 
statarray4 = grpstats(trialTable,{'Subject','HeadTilt','ImageType','ImageTilt'},{'median'},'DataVars','median_T');
statarray5 = grpstats(statarray4,{'Subject','HeadTilt','ImageType'},'mean','DataVars','median_median_T');
ok = repelem(statarray5.mean_median_median_T,3);
statarray4.mean_median_Torsion = ok;

% Head Rotation stuff
% Correct the head roll numbers bc they are inverted. With this correction,
% a positive head roll value means a head roll toward the right shoulder.
% Actually, we will not use this. But just FYI.
saccade.HeadRoll_corrected = -(saccade.HeadRoll);

% Get the median head roll value for each subject for each head tilt
% condition
out = grpstats(saccade,{'Subject','HeadTilt'},{'median'},'DataVars','HeadRoll');

% Put it into the saccade table
HeadTilts = ["Left","Right","Upright"];
for subj = 1:length(SubjList)
    for tilt = 1:length(HeadTilts)
        where = out(find(out.Subject == SubjList(subj) & out.HeadTilt == HeadTilts(tilt)),:);
        saccade.median_HeadRoll(saccade.Subject == SubjList(subj) & saccade.HeadTilt == HeadTilts(tilt)) = table2array(where(1,4)); % the 4 here means just get the fourth col of "out" which is the head roll
    end
end

% Save for future use
MedianHeadTiltRightAmt = -table2array(out(find(out.HeadTilt == 'Right'),4));
MedianHeadTiltLeftAmt = -table2array(out(find(out.HeadTilt == 'Left'),4));
MedianHeadTiltUprightAmt = -table2array(out(find(out.HeadTilt == 'Upright'),4));

% Save to table
tab = table();
tab.Subject = transpose(SubjList);
tab.MedianHeadTiltRightAmt = MedianHeadTiltRightAmt;
tab.MedianHeadTiltLeftAmt = MedianHeadTiltLeftAmt;
tab.MedianHeadTiltUprightAmt = MedianHeadTiltUprightAmt;


%% Implement the circular KDE + circular cross correlation procedure. For fractals + head tilt. 
RightOrLeft = ["Right","Left"];
xcorrKDE_HeadUpHeadRight = [];
xcorrKDE_HeadUpHeadLeft = [];

for aheadtilt = 1:length(RightOrLeft)
    if (DEBUG)
        figure;
        sgtitle(sprintf('Head Upright vs. %s',RightOrLeft(aheadtilt)));
    end
    
    for subj = 1:length(SubjList)
        % Head right vs. Upright for fractals (all images)
        % Get the saccade direction vector that you are interested in
        h1 = saccade.Direction(find(saccade.HeadTilt == RightOrLeft(aheadtilt) & saccade.ImageType == "fractals" & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        h2 = saccade.Direction(find(saccade.HeadTilt == "Upright" & saccade.ImageType == "fractals" & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        
        % Calculate KDE for both distributions
        [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1);  % The .1 is arbitrary (it's the width of the kernel), but seems to work well 
        [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1); 
        
        % Do the circular cross correlation (circular bc we do the
        % repeating)
        hist1Rep = [vfEstimate vfEstimate vfEstimate];
        hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
        [x,lags] = xcorr(hist2Rep, hist1Rep, 450); % 45 degs
        [~,i] = max(x);
        maxlag = lags(i)/10;
        
        % Accumulate these values
        switch RightOrLeft(aheadtilt)
            case "Right"
                xcorrKDE_HeadUpHeadRight = [xcorrKDE_HeadUpHeadRight maxlag];
            case "Left"
                xcorrKDE_HeadUpHeadLeft = [xcorrKDE_HeadUpHeadLeft maxlag];
        end
        
        % Get visualization for each subject
        if (DEBUG)
            subplot(3,14,subj); polarplot([0:.1:360]/180*pi, vfEstimate); hold on; polarplot([0:.1:360]/180*pi, vfEstimate2);
            title(sprintf('Subj %s',SubjList(subj)))
            subplot(3,14,subj+length(SubjList)); plot(vfEstimate); hold on; plot(vfEstimate2);
            subplot(3,14,subj+length(SubjList)*2); plot(lags, x); line([maxlag*10 maxlag*10],get(gca,'ylim')); ...
                text(maxlag, max(x), sprintf('Max @ %.2f deg',maxlag), 'VerticalAlignment','bottom'); ...
                ylim([80 130])
        end
        
    end
end

% Save to table
tab.xcorrKDE_HeadUpHeadRight(:,1) = xcorrKDE_HeadUpHeadRight;
tab.xcorrKDE_HeadUpHeadLeft(:,1) = xcorrKDE_HeadUpHeadLeft;


%% Implement circular KDE and circular cross-corr for fractals in world-referenced coordinates 
clear vfEstimate vfEstimate2 hist1Rep hist2Rep maxlag lags x i

crossCorrFractalsWorldRef_HUHR = [];
crossCorrFractalsWorldRef_HUHL = [];

for aheadtilt = 1:length(RightOrLeft)
    if (DEBUG)
        figure;
        sgtitle(sprintf('Head Upright vs. %s',RightOrLeft(aheadtilt)));
    end
    
    for subj = 1:length(SubjList)
        % Head right vs. Upright for fractals (all images)
        % Get the saccade direction vector that you are interested in
        h1 = saccade.DirAccountingForHeadTilt(find(saccade.HeadTilt == RightOrLeft(aheadtilt) & saccade.ImageType == "fractals" & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        h2 = saccade.Direction(find(saccade.HeadTilt == "Upright" & saccade.ImageType == "fractals" & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        
        % Calculate KDE for both distributions
        [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1)  % 0.01745 radians is 1 deg. The .1 is arbitrary -- maybe change!!
        [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1)  % 0.01745 radians is 1 deg
        
        % Do the circular cross correlation (circular bc we do the
        % repeating)
        hist1Rep = [vfEstimate vfEstimate vfEstimate];
        hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
        [x,lags] = xcorr(hist2Rep, hist1Rep, 450); 
        [~,i] = max(x);
        maxlag = lags(i)/10;
        
        % Accumulate these values
        switch RightOrLeft(aheadtilt)
            case "Right"
                crossCorrFractalsWorldRef_HUHR = [crossCorrFractalsWorldRef_HUHR maxlag];
            case "Left"
                crossCorrFractalsWorldRef_HUHL = [crossCorrFractalsWorldRef_HUHL maxlag];
        end
        
        if (DEBUG)
            subplot(3,14,subj); polarplot([0:.1:360]/180*pi, vfEstimate); hold on; polarplot([0:.1:360]/180*pi, vfEstimate2);
            title(sprintf('Subj %s',SubjList(subj)))
            subplot(3,14,subj+length(SubjList)); plot(vfEstimate); hold on; plot(vfEstimate2);
            subplot(3,14,subj+length(SubjList)*2); plot(lags, x); line([maxlag*10 maxlag*10],get(gca,'ylim')); ...
                text(maxlag, max(x), sprintf('Max @ %.2f deg',maxlag), 'VerticalAlignment','bottom'); ...
                ylim([80 130])
        end

    end
end

% Save to table 
tab.crossCorrFractalsWorldRef_HUHR(:,1) = crossCorrFractalsWorldRef_HUHR;
tab.crossCorrFractalsWorldRef_HUHL(:,1) = crossCorrFractalsWorldRef_HUHL;

%% Quantifying the Offset for Fractals -- bootstrapping error bars for ind subjects 
RightOrLeft = ["Right","Left"];

% Do this for both -30 and 30 image tilts
for aheadtilt = 1:length(RightOrLeft)
    
    % Do the bootstrapping numSampling times
    numSampling = 1000;
    
    % Initialize conglomoration of displacements
    displacements = NaN([length(SubjList),numSampling]);
    
    % Get sacc dirs needed
    for subj = 1:length(SubjList)
        h1 = saccade.Direction(find(saccade.HeadTilt == RightOrLeft(aheadtilt) & saccade.ImageType == "fractals" & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        h2 = saccade.Direction(find(saccade.HeadTilt == "Upright" & saccade.ImageType == "fractals" & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        
        % Get the sampling indicies
        idx_h1 = randi(length(h1),length(h1),numSampling);
        idx_h2 = randi(length(h2),length(h2),numSampling);
        
        % Sample with replacement numSampling times
        for col = 1:numSampling
            newh1 = h1(idx_h1(:,col));
            newh2 = h2(idx_h2(:,col));
            
            % Calculate KDE for both distributions
            [vfEstimate] = circ_ksdensity(newh1,0:.1/180*pi:2*pi,[-pi pi],.1);  % The .1 is arbitrary (it's the width of the kernel), but seems to work well
            [vfEstimate2] = circ_ksdensity(newh2,0:.1/180*pi:2*pi,[-pi pi],.1);
            
            % Do the circular cross correlation (circular bc we do the
            % repeating)
            hist1Rep = [vfEstimate vfEstimate vfEstimate];
            hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
            [x,lags] = xcorr(hist2Rep, hist1Rep, 450); % 45 degs
            [~,i] = max(x);
            maxlag = lags(i)/10;
            
            % Conglomorate
            displacements(subj,col) = maxlag;
        end
    end
    
    
    % Make two tables -- one for 30 and one for -30
    switch aheadtilt
        case 1 % on the first iteration, which is head right
            bootstrapped_right = table();
            bootstrapped_right.Subject = transpose(SubjList);
            bootstrapped_right.SubjNum(:,1) = [1:14];
            bootstrapped_right.Avg(:,1) = 0;
            bootstrapped_right.CI_low(:,1) = 0;
            bootstrapped_right.CI_high(:,1) = 0;
            for subj = 1:length(SubjList)
                bootstrapped_right.Avg(subj,1) = mean(displacements(subj,:));
                bootstrapped_right.CI_low(subj,1) = prctile(displacements(subj,:),2.5);
                bootstrapped_right.CI_high(subj,1) = prctile(displacements(subj,:),97.5);
            end
            
        case 2
            bootstrapped_left = table();
            bootstrapped_left.Subject = transpose(SubjList);
            bootstrapped_left.SubjNum(:,1) = [1:14];
            bootstrapped_left.Avg(:,1) = 0;
            bootstrapped_left.CI_low(:,1) = 0;
            bootstrapped_left.CI_high(:,1) = 0;
            for subj = 1:length(SubjList)
                bootstrapped_left.Avg(subj,1) = mean(displacements(subj,:))
                bootstrapped_left.CI_low(subj,1) = prctile(displacements(subj,:),2.5)
                bootstrapped_left.CI_high(subj,1) = prctile(displacements(subj,:),97.5)
            end
    end
end

% Sort subjects by combined right/left effect 
this = (abs(bootstrapped_left.Avg) + abs(bootstrapped_right.Avg))/2;
[sure,idx] = sortrows(this,'descend')
bootstrapped_right = bootstrapped_right(idx,:);
bootstrapped_left = bootstrapped_left(idx,:);

%% Earth upright scenes at diff head tilt -- calculate KDE distributions + do circular cross correlations
% Do cross correlations for each subject on direction data corrected for head tilt 
RightOrLeft = ["Right","Left"];
crossCorrHeadUprightandRightWithImageEarthUpright = [];
crossCorrHeadUprightandLeftWithImageEarthUpright = [];

for aheadtilt = 1:length(RightOrLeft)
    if (DEBUG)
        figure;
        sgtitle(sprintf('Head Upright vs. %s',RightOrLeft(aheadtilt)));
    end
    
    for subj = 1:length(SubjList)
        
        % Get the saccade direction vector that you are interested in
        switch RightOrLeft(aheadtilt)
            case "Right"
                h1 = saccade.DirAccountingForHeadTilt(find(saccade.HeadTilt == RightOrLeft(aheadtilt) & saccade.ImageType == "scenes" & saccade.ImageTilt == -30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
            case "Left"
                h1 = saccade.DirAccountingForHeadTilt(find(saccade.HeadTilt == RightOrLeft(aheadtilt) & saccade.ImageType == "scenes" & saccade.ImageTilt == 30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        end     
        h2 = saccade.Direction(find(saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.ImageTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        
        % Calculate KDE for both distributions
        [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1)  % The .1 is arbitrary (it's the width of the kernel), but seems to work well 
        [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1)  
        
        % Do the circular cross correlation (circular bc we do the
        % repeating)
        hist1Rep = [vfEstimate vfEstimate vfEstimate];
        hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
        [x,lags] = xcorr(hist2Rep, hist1Rep, 450); % 45 degs
        [~,i] = max(x);
        maxlag = lags(i)/10;
        
        % Accumulate these values
        switch RightOrLeft(aheadtilt)
            case "Right"
                crossCorrHeadUprightandRightWithImageEarthUpright = [crossCorrHeadUprightandRightWithImageEarthUpright maxlag];
            case "Left"
                crossCorrHeadUprightandLeftWithImageEarthUpright = [crossCorrHeadUprightandLeftWithImageEarthUpright maxlag];
        end
        
        % Get visualization for each subject
        if (DEBUG)
            subplot(3,14,subj); polarplot([0:.1:360]/180*pi, vfEstimate); hold on; polarplot([0:.1:360]/180*pi, vfEstimate2);
            title(sprintf('Subj %s',SubjList(subj)))
            subplot(3,14,subj+length(SubjList)); plot(vfEstimate); hold on; plot(vfEstimate2);
            subplot(3,14,subj+length(SubjList)*2); plot(lags, x); line([maxlag*10 maxlag*10],get(gca,'ylim')); ...
                text(maxlag, max(x), sprintf('Max @ %.2f deg',maxlag), 'VerticalAlignment','bottom'); ...
                ylim([80 130])
        end
        
    end
end

% Save to table
tab.crossCorrHeadUprightandRightWithImageEarthUpright(:,1) = crossCorrHeadUprightandRightWithImageEarthUpright;
tab.crossCorrHeadUprightandLeftWithImageEarthUpright(:,1) = crossCorrHeadUprightandLeftWithImageEarthUpright;


%% Earth upright scenes at diff head tilt -- bootstrapping error bars for ind subjects 

% Do this for both -30 and 30 image tilts
for aheadtilt = 1:length(RightOrLeft)
    
    % Do the bootstrapping numSampling times. 1000 in manuscript (it takes
    % a while!!)
    numSampling = 1000;
    
    % Initialize conglomoration of displacements
    displacements = NaN([length(SubjList),numSampling]);
    
    % Get sacc dirs needed
    for subj = 1:length(SubjList)
        switch RightOrLeft(aheadtilt)
            case "Right"
                h1 = saccade.DirAccountingForHeadTilt(find(saccade.HeadTilt == RightOrLeft(aheadtilt) & saccade.ImageType == "scenes" & saccade.Subject == SubjList(subj) & saccade.ImageTilt == -30 & saccade.TrialNumber > 0));
            case "Left"
                h1 = saccade.DirAccountingForHeadTilt(find(saccade.HeadTilt == RightOrLeft(aheadtilt) & saccade.ImageType == "scenes" & saccade.Subject == SubjList(subj) & saccade.ImageTilt == 30 & saccade.TrialNumber > 0));
        end
        h2 = saccade.Direction(find(saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.Subject == SubjList(subj) & saccade.ImageTilt == 0 & saccade.TrialNumber > 0));
        
        % Get the sampling indicies
        idx_h1 = randi(length(h1),length(h1),numSampling);
        idx_h2 = randi(length(h2),length(h2),numSampling);
        
        % Sample with replacement numSampling times
        for col = 1:numSampling
            newh1 = h1(idx_h1(:,col));
            newh2 = h2(idx_h2(:,col));
            
            % Calculate KDE for both distributions
            [vfEstimate] = circ_ksdensity(newh1,0:.1/180*pi:2*pi,[-pi pi],.1);  % The .1 is arbitrary (it's the width of the kernel), but seems to work well
            [vfEstimate2] = circ_ksdensity(newh2,0:.1/180*pi:2*pi,[-pi pi],.1);
            
            % Do the circular cross correlation (circular bc we do the
            % repeating)
            hist1Rep = [vfEstimate vfEstimate vfEstimate];
            hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
            [x,lags] = xcorr(hist2Rep, hist1Rep, 450); % 45 degs
            [~,i] = max(x);
            maxlag = lags(i)/10;
            
            % Conglomorate
            displacements(subj,col) = maxlag;
        end
    end
    
    % Make two tables -- one for 30 and one for -30
    switch aheadtilt
        case 1 % on the first iteration, which is head right
            bootstrapped_30 = table();
            bootstrapped_30.Subject = transpose(SubjList);
            bootstrapped_30.SubjNum(:,1) = [1:14];
            bootstrapped_30.Avg(:,1) = 0;
            bootstrapped_30.CI_low(:,1) = 0;
            bootstrapped_30.CI_high(:,1) = 0;
            for subj = 1:length(SubjList)
                bootstrapped_30.Avg(subj,1) = mean(displacements(subj,:));
                bootstrapped_30.CI_low(subj,1) = prctile(displacements(subj,:),2.5);
                bootstrapped_30.CI_high(subj,1) = prctile(displacements(subj,:),97.5);
            end
            
        case 2
            bootstrapped_neg30 = table();
            bootstrapped_neg30.Subject = transpose(SubjList);
            bootstrapped_neg30.SubjNum(:,1) = [1:14];
            bootstrapped_neg30.Avg(:,1) = 0;
            bootstrapped_neg30.CI_low(:,1) = 0;
            bootstrapped_neg30.CI_high(:,1) = 0;
            for subj = 1:length(SubjList)
                bootstrapped_neg30.Avg(subj,1) = mean(displacements(subj,:));
                bootstrapped_neg30.CI_low(subj,1) = prctile(displacements(subj,:),2.5);
                bootstrapped_neg30.CI_high(subj,1) = prctile(displacements(subj,:),97.5);
            end
    end
end

% Sort subjects by combined right/left effect 
this = (abs(bootstrapped_30.Avg) + abs(bootstrapped_neg30.Avg))/2;
[sure,idx] = sortrows(this,'descend')
bootstrapped_30 = bootstrapped_30(idx,:);
bootstrapped_neg30 = bootstrapped_neg30(idx,:);



%% Image Tilt effect during head upright KDE and cross-correlation
crossCorrHeadUprightImageTilt30and0 = [];
crossCorrHeadUprightImageTiltneg30and0 = [];

ImageTilts = [30,-30];

for animagetilt = 1:length(ImageTilts)
    if (DEBUG)
        figure;
        sgtitle(sprintf('Scene 0 vs. %f',ImageTilts(animagetilt)));
    end
    
    for subj = 1:length(SubjList)
        
        % Get the saccade direction vector that you are interested in
        switch ImageTilts(animagetilt)
            case 30
                h1 = saccade.Direction(find(saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.ImageTilt == 30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
            case -30
                h1 = saccade.Direction(find(saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.ImageTilt == -30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        end     
        h2 = saccade.Direction(find(saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.ImageTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        
        % Calculate KDE for both distributions
        [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1);  % The .1 is arbitrary (it's the width of the kernel), but seems to work well 
        [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1); 
        
        % Do the circular cross correlation (circular bc we do the
        % repeating)
        hist1Rep = [vfEstimate vfEstimate vfEstimate];
        hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
        [x,lags] = xcorr(hist2Rep, hist1Rep, 450); % 45 degs
        [~,i] = max(x);
        maxlag = lags(i)/10;
        
        % Accumulate these values
        switch ImageTilts(animagetilt)
            case 30
                crossCorrHeadUprightImageTilt30and0 = [crossCorrHeadUprightImageTilt30and0 maxlag];
            case -30
                crossCorrHeadUprightImageTiltneg30and0 = [crossCorrHeadUprightImageTiltneg30and0 maxlag];
        end
        
        % Get visualization for each subject
        if (DEBUG)
            subplot(3,14,subj); polarplot([0:.1:360]/180*pi, vfEstimate); hold on; polarplot([0:.1:360]/180*pi, vfEstimate2);
            title(sprintf('Subj %f',SubjList(subj)))
            subplot(3,14,subj+length(SubjList)); plot(vfEstimate); hold on; plot(vfEstimate2);
            subplot(3,14,subj+length(SubjList)*2); plot(lags, x); line([maxlag*10 maxlag*10],get(gca,'ylim')); ...
                text(maxlag, max(x), sprintf('Max @ %.2f deg',maxlag), 'VerticalAlignment','bottom'); ...
                ylim([80 130])
        end
        
    end
end

tab.crossCorrHeadUprightImageTilt30and0(:,1) = crossCorrHeadUprightImageTilt30and0;
tab.crossCorrHeadUprightImageTiltneg30and0(:,1) = crossCorrHeadUprightImageTiltneg30and0;


%% Image Tilt effect during head upright -- bootstrapping error bars for ind subjects 
ImageTilts = [30,-30];

% Do this for both -30 and 30 image tilts
for animagetilt = 1:length(ImageTilts)
    
    % Do the bootstrapping numSampling times
    numSampling = 1000;
    
    % Initialize conglomoration of displacements
    displacements = NaN([length(SubjList),numSampling]);
    
    for subj = 1:length(SubjList)
        h1 = saccade.Direction(find(saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.Subject == SubjList(subj) & saccade.ImageTilt == ImageTilts(animagetilt) & saccade.TrialNumber > 0));
        h2 = saccade.Direction(find(saccade.HeadTilt == "Upright" & saccade.ImageType == "scenes" & saccade.Subject == SubjList(subj) & saccade.ImageTilt == 0 & saccade.TrialNumber > 0));
        
        % Get the sampling indicies
        idx_h1 = randi(length(h1),length(h1),numSampling);
        idx_h2 = randi(length(h2),length(h2),numSampling);
        
        % Sample with replacement numSampling times
        for col = 1:numSampling
            newh1 = h1(idx_h1(:,col));
            newh2 = h2(idx_h2(:,col));
            
            % Calculate KDE for both distributions
            [vfEstimate] = circ_ksdensity(newh1,0:.1/180*pi:2*pi,[-pi pi],.1);  % The .1 is arbitrary (it's the width of the kernel), but seems to work well
            [vfEstimate2] = circ_ksdensity(newh2,0:.1/180*pi:2*pi,[-pi pi],.1);
            
            % Do the circular cross correlation (circular bc we do the
            % repeating)
            hist1Rep = [vfEstimate vfEstimate vfEstimate];
            hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
            [x,lags] = xcorr(hist2Rep, hist1Rep, 450); % 45 degs
            [~,i] = max(x);
            maxlag = lags(i)/10;
            
            % Conglomorate
            displacements(subj,col) = maxlag;
        end
    end
    
    % Make two tables -- one for 30 and one for -30
    switch animagetilt
        case 1 % on the first iteration, which is imagetilt 30
            bootstrapped_30_Im = table();
            bootstrapped_30_Im.Subject = transpose(SubjList);
            bootstrapped_30_Im.SubjNum(:,1) = [1:14];
            bootstrapped_30_Im.Avg(:,1) = 0;
            bootstrapped_30_Im.CI_low(:,1) = 0;
            bootstrapped_30_Im.CI_high(:,1) = 0;
            for subj = 1:length(SubjList)
                bootstrapped_30_Im.Avg(subj,1) = mean(displacements(subj,:));
                bootstrapped_30_Im.CI_low(subj,1) = prctile(displacements(subj,:),2.5);
                bootstrapped_30_Im.CI_high(subj,1) = prctile(displacements(subj,:),97.5);
            end
            
        case 2
            bootstrapped_neg30_Im = table();
            bootstrapped_neg30_Im.Subject = transpose(SubjList);
            bootstrapped_neg30_Im.SubjNum(:,1) = [1:14];
            bootstrapped_neg30_Im.Avg(:,1) = 0;
            bootstrapped_neg30_Im.CI_low(:,1) = 0;
            bootstrapped_neg30_Im.CI_high(:,1) = 0;
            for subj = 1:length(SubjList)
                bootstrapped_neg30_Im.Avg(subj,1) = mean(displacements(subj,:));
                bootstrapped_neg30_Im.CI_low(subj,1) = prctile(displacements(subj,:),2.5);
                bootstrapped_neg30_Im.CI_high(subj,1) = prctile(displacements(subj,:),97.5);
            end
            

    end
end
% Sort subjects by combined right/left effect 
this = (abs(bootstrapped_30_Im.Avg) + abs(bootstrapped_neg30_Im.Avg))/2;
[sure,idx] = sortrows(this,'descend')
bootstrapped_30_Im = bootstrapped_30_Im(idx,:);
bootstrapped_neg30_Im = bootstrapped_neg30_Im(idx,:);


%% Image tilt effect (scenes) for head right and left 
crossCorrHeadRightImageTilt30and0 = [];
crossCorrHeadRightImageTiltneg30and0 = [];

ImageTilts = [30,-30];

for animagetilt = 1:length(ImageTilts)
    if (DEBUG)
        figure;
        sgtitle(sprintf('Scene 0 vs. %f',ImageTilts(animagetilt)));
    end
    
    for subj = 1:length(SubjList)
        
        % Get the saccade direction vector that you are interested in
        switch ImageTilts(animagetilt)
            case 30
                h1 = saccade.Direction(find(saccade.HeadTilt == "Right" & saccade.ImageType == "scenes" & saccade.ImageTilt == 30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
            case -30
                h1 = saccade.Direction(find(saccade.HeadTilt == "Right" & saccade.ImageType == "scenes" & saccade.ImageTilt == -30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        end     
        h2 = saccade.Direction(find(saccade.HeadTilt == "Right" & saccade.ImageType == "scenes" & saccade.ImageTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        
        % Calculate KDE for both distributions
        [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1);  % The .1 is arbitrary (it's the width of the kernel), but seems to work well 
        [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1); 
        
        % Do the circular cross correlation (circular bc we do the
        % repeating)
        hist1Rep = [vfEstimate vfEstimate vfEstimate];
        hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
        [x,lags] = xcorr(hist2Rep, hist1Rep, 450); % 45 degs
        [~,i] = max(x);
        maxlag = lags(i)/10;
        
        % Accumulate these values
        switch ImageTilts(animagetilt)
            case 30
                crossCorrHeadRightImageTilt30and0 = [crossCorrHeadRightImageTilt30and0 maxlag];
            case -30
                crossCorrHeadRightImageTiltneg30and0 = [crossCorrHeadRightImageTiltneg30and0 maxlag];
        end
        
        % Get visualization for each subject
        if (DEBUG)
            subplot(3,14,subj); polarplot([0:.1:360]/180*pi, vfEstimate); hold on; polarplot([0:.1:360]/180*pi, vfEstimate2);
            title(sprintf('Subj %f',SubjList(subj)))
            subplot(3,14,subj+length(SubjList)); plot(vfEstimate); hold on; plot(vfEstimate2);
            subplot(3,14,subj+length(SubjList)*2); plot(lags, x); line([maxlag*10 maxlag*10],get(gca,'ylim')); ...
                text(maxlag, max(x), sprintf('Max @ %.2f deg',maxlag), 'VerticalAlignment','bottom'); ...
                ylim([80 130])
        end
        
    end
end

tab.crossCorrHeadRightImageTilt30and0(:,1) = crossCorrHeadRightImageTilt30and0;
tab.crossCorrHeadRightImageTiltneg30and0(:,1) = crossCorrHeadRightImageTiltneg30and0;


% Do the same as above but for left head tilt... sry i'm a bad programmer
crossCorrHeadLeftImageTilt30and0 = [];
crossCorrHeadLeftImageTiltneg30and0 = [];

ImageTilts = [30,-30];

for animagetilt = 1:length(ImageTilts)
    if (DEBUG)
        figure;
        sgtitle(sprintf('Scene 0 vs. %f',ImageTilts(animagetilt)));
    end
    
    for subj = 1:length(SubjList)
        
        % Get the saccade direction vector that you are interested in
        switch ImageTilts(animagetilt)
            case 30
                h1 = saccade.Direction(find(saccade.HeadTilt == "Left" & saccade.ImageType == "scenes" & saccade.ImageTilt == 30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
            case -30
                h1 = saccade.Direction(find(saccade.HeadTilt == "Left" & saccade.ImageType == "scenes" & saccade.ImageTilt == -30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        end     
        h2 = saccade.Direction(find(saccade.HeadTilt == "Left" & saccade.ImageType == "scenes" & saccade.ImageTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
        
        % Calculate KDE for both distributions
        [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1);  % The .1 is arbitrary (it's the width of the kernel), but seems to work well 
        [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1); 
        
        % Do the circular cross correlation (circular bc we do the
        % repeating)
        hist1Rep = [vfEstimate vfEstimate vfEstimate];
        hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
        [x,lags] = xcorr(hist2Rep, hist1Rep, 450); % 45 degs
        [~,i] = max(x);
        maxlag = lags(i)/10;
        
        % Accumulate these values
        switch ImageTilts(animagetilt)
            case 30
                crossCorrHeadLeftImageTilt30and0 = [crossCorrHeadLeftImageTilt30and0 maxlag];
            case -30
                crossCorrHeadLeftImageTiltneg30and0 = [crossCorrHeadLeftImageTiltneg30and0 maxlag];
        end
        
        % Get visualization for each subject
        if (DEBUG)
            subplot(3,14,subj); polarplot([0:.1:360]/180*pi, vfEstimate); hold on; polarplot([0:.1:360]/180*pi, vfEstimate2);
            title(sprintf('Subj %f',SubjList(subj)))
            subplot(3,14,subj+length(SubjList)); plot(vfEstimate); hold on; plot(vfEstimate2);
            subplot(3,14,subj+length(SubjList)*2); plot(lags, x); line([maxlag*10 maxlag*10],get(gca,'ylim')); ...
                text(maxlag, max(x), sprintf('Max @ %.2f deg',maxlag), 'VerticalAlignment','bottom'); ...
                ylim([80 130])
        end
        
    end
end

tab.crossCorrHeadLeftImageTilt30and0(:,1) = crossCorrHeadLeftImageTilt30and0;
tab.crossCorrHeadLeftImageTiltneg30and0(:,1) = crossCorrHeadLeftImageTiltneg30and0;

%% Save files 
saveTheseTables = ["bootstrapped_30","bootstrapped_neg30","bootstrapped_30_Im",...
    "bootstrapped_neg30_Im","bootstrapped_right","bootstrapped_left","tab"];

if (saveFiles)
    extention='.mat';
    
    for atable = 1:length(saveTheseTables)
        name=saveTheseTables(atable);
        matname = fullfile(saveHere, [char(name) extention]);
        save(matname, saveTheseTables(atable));
    end
    
end













