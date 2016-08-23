% Script to develop blue whale locations 
% This 2nd version tries an alternative way of locating whales (see July 7 notes - option 3) 
% This creates location in structure loc that 
%
% This 3rd version explores a methods of finding the optimal set of travel times based on maximizing the number of
% fitting arrivals times the quality divided by the travel time
% This creates locations in structure loc2
%
% This version superceeds locate3.m

% Adjust to get plots with a nice font size
set(0,'defaultaxesfontsize',18)
set(0,'defaulttextfontsize',18)


%% Parameters
% *** Commented out the next 3 inputs if using monthByMonth_locate3_working ***
% % % File in which to save locations
%  p.fileLocation = 'locations_Jan20_thesh8.mat';           
% % % % Time limits to process
%  p.startDate = datenum('Jan 20 2012');             
%  p.endDate = datenum('Jan 21 2012');

% Input Files
p.fileDetections = 'myCORRECTdetectionsYearA2_2.mat';
p.fileStations = 'usedstations.mat';
eval('load usedstations.mat')

% Plot directory
p.dirPlots = 'December14/likelihood_version3';    % Directory for plots

% Master detections (minimum amplitude, duration and strength (duration*amplitude)
p.masterThreshold = 5; %4;
p.masterDuration = 4; %4;
p.masterStrength = 0;

% Other detections (minimum amplitude, duration and strength (duration*amplitude)
p.minThreshold = 2.2; %1;
p.minDuration = .1; %2;
p.minStrength = 0; %4;

% Maxium range of stations from master stations
p.maxRangeStation = 125;

% Maximum time to look before master detection (the master detection might not be the earliest arrival)
p.maxPreMaster = 0;

% Velocity (km/s) - presently using a constant water velocity
p.velocity = 1.5;

% Minimum number of stations to attempt a location
p.minStation = 4;

% Time and spatial increments for location grid
p.dt = 0.5;
p.dxy = 1;

% Half-width dimensions of grid search in time and space
p.halfWidthOT = 3;
p.halfWidthxy = 100;

% Look up table to get detection timing uncertainty from strength (product
% of amplitude and duration - presently configured to give a constant uncertainty
p.tableStrengthSigma = [0 1; 1e10 1];

% Plot solution (1 - Summed log likelihood; (2 - Also likelihood for individual stations (NOT IMPLEMENTED), 
%                Negative - Print plots to jpeg files
p.plotLikelihood = 0; %-1;

% For the original method, spatial uncertainties greater than this in X or Y lead to a failed status
p.maxSD = 20;

% If the arrival time misfit for a given X-Y-T is > this many multiples of the timing uncertainty is set it to this value
% (mean that detections that fit badly everywhere do not influence preferred solution)
p.maxMisfit = 6;

% For new method, a detection that is fitted to within this many multiples of the timing uncertainty is considered adequately fit
p.okayMisfit = 3;

% For the new method which seeks solutions with strong arrivals near the solution an arrival this close to the call in 
% time is considered optimally fit in time
p.minEarlyFit = 20;



%% Process inputs
% Create single element structure of station locations and spacings
eval(['load ' p.fileStations])
staVec = compileStations(station); 

% Create 3D matrix of times and ranges
[Tpred,xGrid,yGrid,nxGrid,nyGrid,range] = createTpred(staVec,[min(staVec.x) max(staVec.x)], ...
          [min(staVec.y) max(staVec.y)],p.halfWidthxy,p.dxy,p.dxy,p.velocity);

% Create vector of sorted detections
eval(['load ' p.fileDetections]);
det = compileDetections(detections);
% Select based on dates and miminmum criteria
i = det.peakAmp>=p.minThreshold & det.duration>=p.minDuration & ...
    det.strength>=p.minStrength & det.peakTime>=p.startDate & ...
    det.peakTime<=p.endDate;
det.stationName = det.stationName(i);
det.peakAmp = det.peakAmp(i);
det.peakTime = det.peakTime(i);
det.duration = det.duration(i);
det.strength = det.strength(i);
% Track detections that have been used in a solution (so they do not get used in two)
det.inLocation = false(size(det.peakTime));
det.inLocation2 = false(size(det.peakTime));

%% Find detection meeting criterial for master detctions
im = 1;
imVec = find(det.peakAmp>=p.masterThreshold & ...
             det.duration>=p.masterDuration & ...
             det.strength>=p.masterStrength );

%% Loop through potential master calls
il = 0;

% Since this a script, clear previous solutions
clear loc loc2

while im<length(det.peakAmp)
  % Index of next detection meeting Master criteria that is not already in a location
  
%   im = imVec(find(imVec>=im & ~det.inLocation(imVec),1));

  im = imVec(find(imVec>=im,1));
  
  % Get out of while loop if there are no more master detections
  if isempty(im), break, end
  
  %% Detections for current master detection
  % Find detections within feasible time window of master station to create an "event"
  
  
%   j = find((det.peakTime(1:end)-det.peakTime(im))*86400 >= -p.maxPreMaster & ...
%     (det.peakTime(1:end)-det.peakTime(im))*86400 <= p.maxRangeStation/p.velocity & ...
%     ~det.inLocation);
  
 j = find((det.peakTime(1:end)-det.peakTime(im))*86400 >= -p.maxPreMaster & ...
    (det.peakTime(1:end)-det.peakTime(im))*86400 <= p.maxRangeStation/p.velocity);

  % Index of master detection in event
  master = find(j==im);
  
  % Get range of each detection statinon to master station
  range2master = zeros(size(j));
  for i=1:length(j);
    range2master(i)=staVec.spacing(find(strcmpi(staVec.name,det.stationName(j(master)))), ...
                            find(strcmpi(staVec.name,det.stationName(j(i)))));
  end
  
  % Eliminate detections that are infeasible
  %  1. Too far away
  %  2. Repeat for master station
  k = false(size(j));
  for i=1:length(k)
    if i==master || (range2master(i)<p.maxRangeStation && ...
            ~strcmpi(det.stationName(j(master)),det.stationName(j(i))))
      k(i) = true;
    end
  end
  j = j(k);
  range2master = range2master(k);
  
  % Recalculate index of master detection in event 
  master = find(j==im);
  
  %% Only proceed if enough detection stations
  if length(unique(det.stationName(j)))>=p.minStation;
    
    %% Print out detections for current event
    fprintf('   Sta  dTime  Range    Amp Length Strgth       Date       Time\n');
    for i = 1:length(j)
      if i~=master
        fprintf('  %s %6.1f %6.0f %6.2f %6.1f %6.1f    %s\n',char(det.stationName(j(i))), ...
               86400*(det.peakTime(j(i))-det.peakTime(j(master))), ...
               range2master(i),det.peakAmp(j(i)),det.duration(j(i)), ...
               det.strength(j(i)),datestr(det.peakTime(j(i))));
      else
        fprintf('* %s %6.1f %6.1f %6.2f %6.1f %6.1f    %s *\n',char(det.stationName(j(i))), ...
               0,range2master(i),det.peakAmp(j(i)),det.duration(j(i)), ...
               det.strength(j(i)),datestr(det.peakTime(j(i))));
      end
    end
    fprintf('\n');
    
    %% Work on current event
    %Times relative to master detection    % 
    t = (det.peakTime(j)-det.peakTime(j(master)))*86400;
    nt = length(t);
    % Uncertainty of time
    sigma = interp1(p.tableStrengthSigma(:,1),p.tableStrengthSigma(:,2),det.strength(j));
    % Detection quality - currently the same as strength (i.e., amplitude * duration) but could change
%     quality = 2*det.peakAmp(j) + det.duration(j);  
     quality = det.peakAmp(j) .* det.duration(j);  
     
    %Get some indicies
    % Names of stations corresponding to detections in this event
    name = det.stationName(j);
    %Indicies of detections within stationVec
    istaVec = zeros(nt,1);
    for i=1:nt
      istaVec(i) = find(strcmp(staVec.name,name(i)));
    end
    %Indicies within staVec for a subset of stations that have detection in this event
    istaLoc = sort(unique(istaVec));
    nstaLoc = length(istaLoc);
    %Indicies of detections within subset of stations that have detection in this event
    iSub = zeros(size(t));
    for i = 1:length(istaLoc);
      iSub = iSub + i*(istaVec==istaLoc(i));
    end
    
    % Detection quality measurements normalized to sum of quality for that OBS
    % Not presently used but the idea is that one would prefer to incorporate strong detections
    for i=1:nt
      qualityNorm(i) = quality(i) / sum(quality(iSub==iSub(i)));
    end

    %Subgrid for locations (only look within p.halfWidthxy of the master detection station
    xmin = floor((staVec.x(istaVec(master))-p.halfWidthxy)/p.dxy)*p.dxy;
    xmax = ceil((staVec.x(istaVec(master))+p.halfWidthxy)/p.dxy)*p.dxy;
    ixGrid = find(xGrid>=xmin & xGrid<=xmax);
    xGridSub = xGrid(ixGrid);
    nxGridSub = length(ixGrid);
    ymin = floor((staVec.y(istaVec(master))-p.halfWidthxy)/p.dxy)*p.dxy;
    ymax = ceil((staVec.y(istaVec(master))+p.halfWidthxy)/p.dxy)*p.dxy;
    iyGrid = find(yGrid>=ymin & yGrid<=ymax);
    yGridSub = yGrid(iyGrid);
    nyGridSub = length(iyGrid);

    %Travel times & Range for subgrid
    TpredSub = Tpred(ixGrid,iyGrid,istaLoc);
%     rangeSub = range(ixGrid,iyGrid,istaLoc);
 
    % Matricies to hold likelyhood, travel time misfit and logical "fit" status for each detection
    pBuild = zeros(nxGridSub,nyGridSub,nt);
    misfit = zeros(nxGridSub,nyGridSub,nt);
    fit = false(nxGridSub,nyGridSub,nt);

    %Origin time grid assuming that master call is correct (origin time assuming master call fit perfectly)
    ot0 = t(master) - TpredSub(:,:,iSub(master));

    %Origin times to consider (relative to master call prediction)
    ot = -p.halfWidthOT:p.dt:p.halfWidthOT;
    not = length(ot);

    % Set up matricies for solution space (probability, measure of fitting quality arrivals early, number of fitting)
    pGrid = zeros(nxGridSub,nyGridSub,not);
    qualityEarlyFit = zeros(nxGridSub,nyGridSub,not);
    nFit = zeros(nxGridSub,nyGridSub,not);

    %Compute locations
    % Cycle through each origin time
    for iot = 1:not
      % Cycle through each detection
      for it = 1:nt
        % Misfit in units of travel time uncertainty
        misfit(:,:,it) = abs(t(it)-TpredSub(:,:,iSub(it))-ot0-ot(iot)) / sigma(it);
        misfit(:,:,it) = min(misfit(:,:,it),p.maxMisfit);
        % Convert misfits to something proportional to the log of the probability
        pBuild(:,:,it) = -1/2 * misfit(:,:,it).^2;
        % Compute a quanity that is maximized for high quality arrivals that are nearby 
        % (but limit nearby to a minimum travel time)
        qualityEarlyFit(:,:,iot) = qualityEarlyFit(:,:,iot) + (misfit(:,:,it)<=p.okayMisfit) * quality(it) ./ ...
                                    max(abs(t(it)-ot0-ot(iot)),p.minEarlyFit);
      end
      
      % Number of fitting times for this origin time
      nFit(:,:,iot) = sum(misfit<p.okayMisfit,3);
      
      % Cycle through stations taking the highest probability arrival
      for is = 1:nstaLoc
        i = find(iSub==is);
        if length(i)==1;
          pGrid(:,:,iot) = pGrid(:,:,iot) + pBuild(:,:,i);
        else
          pGrid(:,:,iot) = pGrid(:,:,iot) + max(pBuild(:,:,i),[],3);
        end
      end
    end

    % Do not care about origin time so sum the probability for all origin times.
    pGridSum = sum(exp(pGrid),3);
    % Normalize probability so it sums to 1.
    pGridSumNorm = pGridSum/sum(sum(pGridSum));
    
    % Need at least p.minStation fitting arrivals so set qualityEarlyFit to zero where there ar not
    qualityEarlyFit = qualityEarlyFit .* (nFit>=p.minStation);
    
    % Do not proceed to alternative loc2 solution if there is no X-Y-OT value with enough fitting detecitons
    if ~any(qualityEarlyFit(:));
      noSolution = true;
      
    else
      noSolution = false;

      %% Now doing the work for an alternative loc2 solution
      % Find indicies of maximum of qualityEarlyFit
      [temp1,i3] = max(qualityEarlyFit,[],3);
      [temp2,i2] = max(temp1,[],2);
      [~,i1] = max(temp2);
      i2 = i2(i1);
      i3 = i3(i1,i2);

      % Determine which times to use based on which fit for the largest value of qualityEarlyFit
      use = zeros(nt,1);
      for it = 1:nt
        use(it) = abs(t(it)-TpredSub(i1,i2,iSub(it))-ot0(i1,i2)-ot(i3)) / sigma(it);
      end
      % Only want to use one time per station 
      % (will very occasionally be two if detections are close enough for 2 to fit)
      for i=1:max(iSub)
        repeat = find(iSub==i & use<=p.okayMisfit);
        if length(repeat)>1
          [~,best] = max(det.strength(j(repeat)));
          repeat = repeat([1:best-1 best+1:end]);
          use(repeat) = inf;
        end
      end      
      use = use<=p.okayMisfit;
      usestatus=find(use);
      % Index of fitting detections
      j2 = j(use);
      % Get the master detection index amongst fitting detections
      master2 = find(j2==im);
      
      % Get subset of times to use
      t2 = (det.peakTime(j2)-det.peakTime(j2(master2)))*86400;
      nt2 = length(t2);
      quality2 = det.peakAmp(j2) .* det.duration(j2);  % Measure of detection quality
      sigma2 = interp1(p.tableStrengthSigma(:,1),p.tableStrengthSigma(:,2),det.strength(j2));
      
      for qualind=1:length(quality2);
          if quality2(qualind) <= 3;
          sigma2(qualind)=4;    
          elseif quality2(qualind) <= 5;
          sigma2(qualind)=3;
          elseif quality2(qualind) <= 8;
              sigma2(qualind)=2;
          end
          
      end
      
      % Get some indicies
      % Names of stations corresponding to detections in this event
      
      name2 = det.stationName(j2);
      thresh2=det.peakAmp(j2);
      [sortthresh,indsthresh]=sort(thresh2,'descend');
      % Indicies of detections within stationVec
      
      istaVec2 = zeros(nt2,1);
      for i=1:nt2
          istaVec2(i) = find(strcmp(staVec.name,name2(i)));
      end
      %Indicies within staVec for a subset of stations that have detection in this event
      istaLoc2 = sort(unique(istaVec2));
      nstaLoc2 = length(istaLoc2);
      %Indicies of detections within subset of stations that have detection in this event
      iSub2 = zeros(size(t2));
      for i = 1:length(istaLoc2);
        iSub2 = iSub2 + i*(istaVec2==istaLoc2(i));
      end
      % Double check that no station is repeated (probably could delete this)
      if ~any(diff(sort(iSub2)))
        disp('Repeated station in final set of times')
        keyboard
      end

      %Travel times & Range for subgrid
      TpredSub2 = Tpred(ixGrid,iyGrid,istaLoc2);

      % Matricies to hold likelyhood and travel time misfit
      pBuild2 = zeros(nxGridSub,nyGridSub,nt2);
      misfit2 = zeros(nxGridSub,nyGridSub,nt2);

      % Probability matrix for solution space
      pGrid2 = zeros(nxGridSub,nyGridSub,not);

      %Compute likelihood grid
      for iot = 1:not
        for it = 1:nt2
          misfit2(:,:,it) = abs(t2(it)-TpredSub2(:,:,iSub2(it))-ot0-ot(iot)) / sigma2(it);
          misfit2(:,:,it) = min(misfit2(:,:,it),p.maxMisfit);
          pBuild2(:,:,it) = -1/2 * misfit2(:,:,it).^2;
        end

        for is = 1:nstaLoc2
          i = find(iSub2==is);
          if length(i)==1;
            pGrid2(:,:,iot) = pGrid2(:,:,iot) + pBuild2(:,:,i);
          else
            % Should never get here since it means that there are two detections for a station in the solution
            disp('Problems')
            keyboard
            pGrid2(:,:,iot) = pGrid2(:,:,iot) + max(pBuild2(:,:,i),[],3);
          end
        end
      end

      % Likelihood (summed and normalized)
      pGrid2Sum = sum(exp(pGrid2),3);
      pGrid2SumNorm = pGrid2Sum/sum(sum(pGrid2Sum));
    end
    
    %Find maximum of qualityEarlyFit - already done above
    junk = squeeze(max(qualityEarlyFit,[],3));
    [dum,imax] = max(junk,[],1);
    [QEmax,jmax]= max(dum);
    imax = imax(jmax);

    %Plot likelihood surface(s) and save jpegs if requested
    if abs(p.plotLikelihood)>=1
      for i=1:4-noSolution
        figure(100+i);
        clf
        orient tall
        if i==1
          imagesc(xGridSub,yGridSub,log(pGridSum)');
          hold on
          plot(xGridSub(imax),yGridSub(jmax),'+w','markersize',18)
          titlestring = 'pGridSum';
        elseif i==2
          imagesc(xGridSub,yGridSub,log(pGrid2Sum)');
          hold on
          plot(xGridSub(imax),yGridSub(jmax),'+w','markersize',18)
          titlestring = 'pGridSum2';
        elseif i==3
          imagesc(xGridSub,yGridSub,(squeeze(max(nFit,[],3)))'.*(squeeze(max(nFit,[],3))'>=p.minStation));
          hold on
          plot(xGridSub(imax),yGridSub(jmax),'+w','markersize',18)
          titlestring = 'Number that Fit';
        elseif i==4
          imagesc(xGridSub,yGridSub,squeeze(max(qualityEarlyFit,[],3))');
          hold on
          plot(xGridSub(imax),yGridSub(jmax),'+w','markersize',18)
          titlestring = 'Quality-Early Fit';
        end
        hold on
        colormap jet
        axis equal
        colorbar('southoutside')
        set(gca,'ydir','normal');
        lim = [xlim ylim];
        for k=2:length(staVec.x)
          if staVec.x(k)>lim(1) && staVec.x(k)<lim(2) && ...
             staVec.y(k)>lim(3) && staVec.y(k)<lim(4)
            plot(staVec.x(k),staVec.y(k),'ok','markerfacecolor','w')
            if any(istaLoc==k)
              text(staVec.x(k),staVec.y(k),staVec.name(k),'color','w')
            end
          end
        end
        title([titlestring '  ' datestr(det.peakTime(j(master)))])
        xlabel('X, km')
        ylabel('Y, km')
%         xlim([450 800])
        if p.plotLikelihood<0
          if i==1
            eval(['print -dpng ' p.dirPlots '/' datestr(det.peakTime(j(master)),'yyyymmddTHHMMSS') '_likelihood.png'])
          elseif i==2
            eval(['print -dpng ' p.dirPlots '/' datestr(det.peakTime(j(master)),'yyyymmddTHHMMSS') '_likelihood2.png'])
          elseif i==3
            eval(['print -dpng ' p.dirPlots '/' datestr(det.peakTime(j(master)),'yyyymmddTHHMMSS') '_nFit.png'])
          elseif i==4
            eval(['print -dpng ' p.dirPlots '/' datestr(det.peakTime(j(master)),'yyyymmddTHHMMSS') '_qualityEarlyFit.png'])
          end            
        end
      end
    end
    
    %% Process solution
    % loc structure is the original solution and loc2 is the new qualityEarlyFit one (alternate)

    % Maximum likelihood gives location
    il = il+1;
    [pGridMax,i3] = max(pGrid,[],3);
    [temp,i2] = max(pGridMax,[],2);
    [pMax,i1] = max(temp);
    i2 = i2(i1);
    i3 = i3(i1,i2);
    loc(il).status = 0;
    loc(il).x = xGridSub(i1);
    loc(il).y = yGridSub(i2);
    loc(il).ot = det.peakTime(j(master))+(ot0(i1,i2)+ot(i3))/86400;
    if ~noSolution
      [pGrid2Max,i3Alt] = max(pGrid2,[],3);
      [temp,i2Alt] = max(pGrid2Max,[],2);
      [p2Max,i1Alt] = max(temp);
      i2Alt = i2Alt(i1Alt);
      i3Alt = i3Alt(i1Alt,i2Alt);
      loc2(il).status = 0;
      loc2(il).x = xGridSub(i1Alt);
      loc2(il).y = yGridSub(i2Alt);
      loc2(il).ot = det.peakTime(j2(master2))+(ot0(i1Alt,i2Alt)+ot(i3Alt))/86400;
      loc2(il).QEmax = QEmax;
    else
      loc2(il).status = 10000;
      loc2(il).QEmax = QEmax;
    end
    
    % Standard deviation of location (normalized probability function) in X and Y
    temp = sum(pGridSumNorm,2);
    loc(il).xmean = sum(temp.*xGridSub(:));
    loc(il).xsd = sqrt(sum(temp.*(xGridSub(:)-loc(il).xmean).^2));
    temp = sum(pGridSumNorm,1);
    loc(il).ymean = sum(temp.*yGridSub);
    loc(il).ysd = sqrt(sum(temp.*(yGridSub-loc(il).ymean).^2));
    if loc(il).xsd>p.maxSD || loc(il).ysd>p.maxSD
      loc(il).status = 1;
    end
    if ~noSolution
      temp = sum(pGrid2SumNorm,2);
      loc2(il).xmean = sum(temp.*xGridSub(:));
      loc2(il).xsd = sqrt(sum(temp.*(xGridSub(:)-loc(il).xmean).^2));
      temp = sum(pGrid2SumNorm,1);
      loc2(il).ymean = sum(temp.*yGridSub);
      loc2(il).ysd = sqrt(sum(temp.*(yGridSub-loc(il).ymean).^2));
      if loc2(il).xsd>p.maxSD || loc(il).ysd>p.maxSD
        loc2(il).status = loc2(il).status + 1;
      end
    end
    
    % Confidence 50% and 95% confidence limits based on contouring normalized probaility
    temp = sort(pGridSumNorm(:),'descend');
% Following sometimes fails since some grid values are essential 0
%     p95 = interp1(cumsum(temp),temp,0.95);
%     p50 = interp1(cumsum(temp),temp,0.5);
    p95 = temp(find(cumsum(temp)>=0.95,1,'first'));    
    p50 = temp(find(cumsum(temp)>=0.5,1,'first'));    
    if ~noSolution
      temp = sort(pGrid2SumNorm(:),'descend');
      p95_2 = temp(find(cumsum(temp)>=0.95,1,'first'));    
      p50_2 = temp(find(cumsum(temp)>=0.5,1,'first'));    
    end
    % Process contours
    loc(il).xmin50 = Inf;
    loc(il).xmax50 = -Inf;
    loc(il).ymin50 = Inf;
    loc(il).ymax50 = -Inf;
    loc(il).xmin95 = Inf;
    loc(il).xmax95 = -Inf;
    loc(il).ymin95 = Inf;
    loc(il).ymax95 = -Inf;
    cs = contourc(xGridSub,yGridSub,pGridSumNorm',[p50 p95]);
    n95 = 0;
    i = 1; 
    while i<size(cs,2);
      xc = cs(1,i+1:i+cs(2,i));
      yc = cs(2,i+1:i+cs(2,i));
      if cs(1,i)==p95
        n95 = n95+1;
        loc(il).xmin95 = min(loc(il).xmin95,min(xc));
        loc(il).xmax95 = max(loc(il).xmax95,max(xc));
        loc(il).ymin95 = min(loc(il).ymin95,min(yc));
        loc(il).ymax95 = max(loc(il).ymax95,max(yc));
      else
        loc(il).xmin50 = min(loc(il).xmin50,min(xc));
        loc(il).xmax50 = max(loc(il).xmax50,max(xc));
        loc(il).ymin50 = min(loc(il).ymin50,min(yc));
        loc(il).ymax50 = max(loc(il).ymax50,max(yc));
      end
      i = i+cs(2,i)+1;
    end
    loc(il).dx95 = loc(il).xmax95-loc(il).xmin95;
    loc(il).dy95 = loc(il).ymax95-loc(il).ymin95;
    % Multiple minima within 95% confidence limits
    if n95>1
      loc(il).status = loc(il).status + 100;
    end
    % 95% condidence limits extend to edge of grid
    if loc(il).xmin95<=min(xGridSub) || loc(il).xmax95>=max(xGridSub) || ...
       loc(il).ymin95<=min(yGridSub) || loc(il).ymax95>=max(yGridSub)
      loc(il).status = loc(il).status + 10;
    end    
    % Process contours for alternate
    if ~noSolution
      loc2(il).xmin50 = Inf;
      loc2(il).xmax50 = -Inf;
      loc2(il).ymin50 = Inf;
      loc2(il).ymax50 = -Inf;
      loc2(il).xmin95 = Inf;
      loc2(il).xmax95 = -Inf;
      loc2(il).ymin95 = Inf;
      loc2(il).ymax95 = -Inf;
      cs = contourc(xGridSub,yGridSub,pGrid2SumNorm',[p50 p95]);
      n95 = 0;
      i = 1; 
      while i<size(cs,2);
        xc = cs(1,i+1:i+cs(2,i));
        yc = cs(2,i+1:i+cs(2,i));
        if cs(1,i)==p95
          n95 = n95+1;
          loc2(il).xmin95 = min(loc2(il).xmin95,min(xc));
          loc2(il).xmax95 = max(loc2(il).xmax95,max(xc));
          loc2(il).ymin95 = min(loc2(il).ymin95,min(yc));
          loc2(il).ymax95 = max(loc2(il).ymax95,max(yc));
        else
          loc2(il).xmin50 = min(loc2(il).xmin50,min(xc));
          loc2(il).xmax50 = max(loc2(il).xmax50,max(xc));
          loc2(il).ymin50 = min(loc2(il).ymin50,min(yc));
          loc2(il).ymax50 = max(loc2(il).ymax50,max(yc));
        end
        i = i+cs(2,i)+1;
      end
      loc2(il).dx95 = loc2(il).xmax95-loc2(il).xmin95;
      loc2(il).dy95 = loc2(il).ymax95-loc2(il).ymin95;
      % Multiple minima within 95% confidence limits
      if n95>1
        loc2(il).status = loc2(il).status + 100;
      end
      % 95% condidence limits extend to edge of grid
      if loc2(il).xmin95<=min(xGridSub) || loc2(il).xmax95>=max(xGridSub) || ...
         loc2(il).ymin95<=min(yGridSub) || loc2(il).ymax95>=max(yGridSub)
        loc2(il).status = loc2(il).status + 10;
      end    
    end
    
    % Predicted times & misfit
%     tMod = det.peakTime(j(master))+(ot0(i1,i2)+ot(i3))/86400+squeeze(TpredSub(i1,i2,iSub));
    tMod = ot0(i1,i2)+ot(i3)+squeeze(TpredSub(i1,i2,iSub));
    tMisfit = (t - tMod);
    misfit = tMisfit./sigma;
    ifit = zeros(length(t),1);
    for i=1:length(istaLoc);
      index = find(iSub==i);
      [temp,k] = min(abs(misfit(index)));
      if abs(temp)<=p.okayMisfit
        ifit(index(k))=i;
      end
    end
    loc(il).nsta = sum(~~ifit);
    if loc(il).nsta<p.minStation
      loc(il).status = loc(il).status+1000;
    end
    loc(il).rmsNorm = sqrt(sum(misfit(~~ifit).^2)/loc(il).nsta);
    loc(il).rms = sqrt(sum(tMisfit(~~ifit).^2)/loc(il).nsta);
    if ~loc(il).status
      det.inLocation(j(~~ifit)) = true;
    end
    if ~noSolution
      tMod2 = ot0(i1Alt,i2Alt)+ot(i3Alt)+squeeze(TpredSub2(i1Alt,i2Alt,iSub2));
      tMisfit2 = (t2 - tMod2);
      misfit2 = tMisfit2./sigma2;
      ifit2 = zeros(length(t2),1);
      for i=1:length(istaLoc2);
        index = find(iSub2==i);
        [temp,k] = min(abs(misfit2(index)));
        if abs(temp)<=p.okayMisfit
          ifit2(index(k))=i;
        end
      end
      ifit2inds=find(ifit2);
      loc2(il).nsta = sum(~~ifit2);
      if loc2(il).nsta<p.minStation
        loc2(il).status = loc2(il).status+1000;
      end
      loc2(il).rmsNorm = sqrt(sum(misfit2(~~ifit2).^2)/loc2(il).nsta);
      loc2(il).rms = sqrt(sum(tMisfit2(~~ifit2).^2)/loc2(il).nsta);
      if ~loc2(il).status
        det.inLocation2(j) = false   
        det.inLocation2(j(usestatus(ifit2inds))) = true;
      end
    end
    
    % Characteristics of the detections in location
    loc(il).medPeakAmp = median(det.peakAmp(j(~~ifit)));
    loc(il).medDuration = median(det.duration(j(~~ifit)));
    loc(il).medStrength = median(det.strength(j(~~ifit)));
    loc(il).meanPeakAmp = mean(det.peakAmp(j(~~ifit)));
    loc(il).meanDuration = mean(det.duration(j(~~ifit)));
    loc(il).meanStrength = mean(det.strength(j(~~ifit)));
    loc(il).maxPeakAmp = max(det.peakAmp(j(~~ifit)));
    loc(il).maxDuration = max(det.duration(j(~~ifit)));
    loc(il).maxStrength = max(det.strength(j(~~ifit)));   
    if ~noSolution
      loc2(il).medPeakAmp = median(det.peakAmp(j(usestatus(ifit2inds))));
      loc2(il).medDuration = median(det.duration(j(usestatus(ifit2inds))));
      loc2(il).medStrength = median(det.strength(j(usestatus(ifit2inds))));
      loc2(il).meanPeakAmp = mean(det.peakAmp(j(usestatus(ifit2inds))));
      loc2(il).meanDuration = mean(det.duration(j(usestatus(ifit2inds))));
      loc2(il).meanStrength = mean(det.strength(j(usestatus(ifit2inds))));
      loc2(il).maxPeakAmp = max(det.peakAmp(j(usestatus(ifit2inds))));
      loc2(il).maxDuration = max(det.duration(j(usestatus(ifit2inds))));
      loc2(il).maxStrength = max(det.strength(j(usestatus(ifit2inds))));
    end
    
    % Information about detections
    loc(il).det.stationName = det.stationName(j);
    loc(il).det.peakTime = det.peakTime(j);
    loc(il).det.peakTimePred = loc(il).det.peakTime(master)+tMod/86400;
    loc(il).det.peakAmp = det.peakAmp(j);
    loc(il).det.duration = det.duration(j);
    loc(il).det.strength = det.strength(j);
    loc(il).det.inLocation = det.inLocation(j);
    loc(il).det.master = master;
    loc(il).det.t = t;
    loc(il).det.tPred = tMod;
    loc(il).det.sigma = sigma;
    loc(il).det.tMisfit = tMisfit;
    loc(il).det.misfit = misfit;
    if ~noSolution
      loc2(il).det.stationName = det.stationName(j);
      loc2(il).det.peakTime = det.peakTime(j);
      loc2(il).det.peakTimePred = loc2(il).det.peakTime(master)+tMod2/86400;
      loc2(il).det.peakAmp = det.peakAmp(j);
      loc2(il).det.duration = det.duration(j);
      loc2(il).det.strength = det.strength(j);
      loc2(il).det.inLocation = det.inLocation2(j);
      loc2(il).det.master = master;
      loc2(il).det.t = t;
      loc2(il).det.tPred = tMod2;
      loc2(il).det.sigma = sigma;
      loc2(il).det.tMisfit2 = tMisfit2;
      loc2(il).det.misfit2 = misfit2;
      loc2(il).det.stationName(loc2(il).det.inLocation)
      
%       det.inLocation2(j)=false;
%       det.inLocation(j)=false
      
    end      
  end
  
  % Increment for next master call
  im = im+1;
end

%% Save output
% locAttempt are all the original solutions including those that fail
if exist('loc') == 1
locAttempt = loc;
% loc is just those that work
loc = loc(~[loc.status]);
% locBest are the best solutions - should make the conditions parameters
% meanStrength>10 is about equivalent to meanDuration>5
locBest =loc([loc.nsta]>=6 & [loc.meanDuration]>5);

% loc2Attempt are all the alternate solutions including those that fail
loc2Attempt = loc2;
% loc2 are alternate solutions that work 
loc2 = loc2(~[loc2.status]);

eval(['save ' ['/burr1/CI/locations/' p.fileLocation] ' loc locBest locAttempt loc2 loc2Attempt p staVec']);
end
%% Plotting solutions - not used anymore
% %% Plot Locations
% figure(201); clf
% plot([staVec.x(2:10)],[staVec.y(2:10)],'sk','markerfacecolor','k','markersize',9)
% hold on
% col = LinColor([loc.ot],[min([loc.ot]) max([loc.ot])],'jet');
% for i=1:length(loc)
%   plot(loc(i).x,loc(i).y,'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:));
% end
% axis equal
% caxis([min([loc.ot]) max([loc.ot])]-floor(min([loc.ot])))
% colorbar
% title('Locations')
% 
% figure(202); clf
% plot([staVec.x(2:10)],[staVec.y(2:10)],'sk','markerfacecolor','k','markersize',9)
% hold on
% col = LinColor([locBest.ot],[min([locBest.ot]) max([locBest.ot])],'jet');
% for i=1:length(locBest)
%   plot(locBest(i).x,locBest(i).y,'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:));
% end
% axis equal
% caxis([min([locBest.ot]) max([locBest.ot])]-floor(min([locBest.ot])))
% colorbar
% title('Best Locations')
% 
% figure(203); clf
% plot([staVec.x(2:10)],[staVec.y(2:10)],'sk','markerfacecolor','k','markersize',9)
% hold on
% col = LinColor([loc2.ot],[min([loc2.ot]) max([loc2.ot])],'jet');
% for i=1:length(loc2)
%   plot(loc2(i).x,loc2(i).y,'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:));
% end
% axis equal
% caxis([min([loc2.ot]) max([loc2.ot])]-floor(min([loc2.ot])))
% colorbar
% title('Alt Locations')
% figure(204); clf
% plot([staVec.x(2:10)],[staVec.y(2:10)],'sk','markerfacecolor','k','markersize',9)
% hold on
% col = LinColor([loc2Best.ot],[min([loc2Best.ot]) max([loc2Best.ot])],'jet');
% for i=1:length(loc2Best)
%   plot(loc2Best(i).x,loc2Best(i).y,'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:));
% end
% axis equal
% caxis([min([loc2Best.ot]) max([loc2Best.ot])]-floor(min([loc2Best.ot])))
% colorbar
% title('Alt Best Locations')
