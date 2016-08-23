% Script to run locations 1 month at a time using locate3_working and create daily plots
% Some plot limits are hardwired
% *** Note that for this to work 3 parameters need to be commented out in locate3_working ***

%% Parameters
% Criteria to group calls in a track
drMax = 10;   %Maximum allowed spacing in kilometers
dtMax = 1/24; %Maximum spacing in time (days)
n4track = 4;  %Number of other calls nearby in space and time to consider a call to be in a track

% Number of calls in a track for it to be signitifant (which means it gets plotted)
nSignificant = 5;

% directory for plots
direc = '/burr1/CI/locations/monthly_locations3';

% Files in wich to save solutions for each interval
file = {'monthly_locations3/mylocation4_Sep2011.mat'; ...
        'monthly_locations3/mylocation4_Oct2011.mat'; ...
        'monthly_locations3/mylocation4_Nov2011.mat'; ...
        'monthly_locations3/mylocation4_Dec2011.mat'; ...
        'monthly_locations3/mylocation4_Jan2012.mat'; ...
        'monthly_locations3/mylocation4_Feb2012.mat'; ...
        'monthly_locations3/mylocation4_Mar2012.mat'; ...
        'monthly_locations3/mylocation4_Apr2012.mat'; ...
        'monthly_locations3/mylocation4_May2012.mat'; ...
        'monthly_locations3/mylocation4_Jun2012.mat'; ...
        'monthly_locations3/mylocation4_Jul2012.mat'};
      
% Time intervals to process
month = { 'September' 'October' 'November' 'December' 'January' 'February' 'March' 'April' 'May' 'June' 'July'};
year = { '2011' '2011' '2011' '2011' '2012' '2012' '2012' '2012' '2012' '2012' '2012'};
startDate = [ datenum('Sep 1 2011') datenum('Oct 1 2011') datenum('Nov 1 2011') datenum('Dec 1 2011') datenum('Jan 1 2012')...
  datenum('Feb 1 2012') datenum('Mar 1 2012') datenum('Apr 1 2012') datenum('May 1 2012') datenum('Jun 1 2012') datenum('Jul 1 2012')];
endDate = [datenum('Oct 1 2011') datenum('Nov 1 2011') datenum('Dec 1 2011') datenum('Jan 1 2012')...
  datenum('Feb 1 2012') datenum('Mar 1 2012') datenum('Apr 1 2012') datenum('May 1 2012')...
  datenum('Jun 1 2012') datenum('Jul 1 2012') datenum('Aug 1 2012')];

%% Loop over intervals to process
for loop = 2:length(file)
  
  % Locate the calls 
  p.fileLocation = file{loop};
  p.startDate = startDate(loop);
  p.endDate = endDate(loop);
  %locate3_working
  
  % Load locations (unecessary at present since locate3_working is a script)
eval(['load ' ['/burr1/CI/locations/' file{loop}]]);
  
  % Filter loc2 locations for tracks
  if exist('loc2') == 1
  loc2Filt = filter4track(loc2,drMax,dtMax,n4track);
  ot = [loc2.ot]-datenum([month{loop} ' 1, ' year{loop}])+1;
  x = [loc2.x];
  y = [loc2.y];
  otFilt = [loc2Filt.ot]-datenum([month{loop} ' 1, ' year{loop}])+1;
  xFilt = [loc2Filt.x];
  yFilt = [loc2Filt.y];
  load('CascadiaBathymetry.mat');
  %% Make daily plots
  for day = 1:31
    
    % Plot loc2 locations
    j = find(ot>=day & ot<day+1);
    if ~isempty(j) 
      figure(201); clf
      contour(Cascadia.xgrid,Cascadia.ygrid,Cascadia.elevgrid);
      hold on
      plot([staVec.x(1:end)],[staVec.y(1:end)],'sk','markerfacecolor','k','markersize',9)
      for i=[40 42 36 31 23 4 12 16 13]
        text(staVec.x(i),staVec.y(i),staVec.name(i),'fontsize',12)
      end
      
      if length(j)>1
        col = LinColor((ot(j)-floor(min(ot(j(1)))))*24,[0 24],'jet');
      else
        col = [0 0 1];
      end
      for i=1:length(j)
        plot(x(j(i)),y(j(i)),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:));
      end
      axis equal
      if length(j)>1
        caxis([0 24]);
        colorbar
      end
      % Hardwire plot limits so each day has the same scale
      xlim([min(staVec.x)-10 max(staVec.x)+50])
      ylim([min(staVec.y)-10 max(staVec.y)+10])
      eval(['title('' Locate3 Working - ' month{loop} ' ' int2str(day) ', ', year{loop} ''',''fontsize'',12)'])
      ylabel('Y, km')
      xlabel('X, km')
      eval(['print -dpng ' direc '/my4_locate3_working' datestr(startDate(loop)+day-1,'yyyy-mm-dd') '.png'])
      %savefig([direc '/my_locate3_working' datestr(startDate(loop)+day-1,'yyyy-mm-dd')])

%       eval(['print -depsc ' direc '\locate3_working_' month{loop} int2str(day) '.eps'])
    end
    
    % plot loc2 locations after filtering for tracks
    j = find(otFilt>=day & otFilt<day+1);
    if ~isempty(j) && length(j)>nSignificant
      figure(202); clf
      contour(Cascadia.xgrid,Cascadia.ygrid,Cascadia.elevgrid);
      hold on
      plot([staVec.x(1:end)],[staVec.y(1:end)],'sk','markerfacecolor','k','markersize',9)
      for i=[40 42 36 31 23 4 12 16 13]
        text(staVec.x(i),staVec.y(i),staVec.name(i),'fontsize',12)
      end
      hold on
      if length(j)>1
        col = LinColor((otFilt(j)-floor(min(otFilt(j(1)))))*24,[0 24],'jet');
      else
        col = [0 0 1];
      end
      for i=1:length(j)
        plot(xFilt(j(i)),yFilt(j(i)),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:));
      end
      axis equal
      if length(j)>1
        caxis([0 24]);
        colorbar
      end
      % Hardwire plot limits so each day has the same scale
      xlim([min(staVec.x)-10 max(staVec.x)+50])
      ylim([min(staVec.y)-10 max(staVec.y)+10])
      eval(['title('' Locate3 Working Filtered - ' month{loop} ' ' int2str(day) ', ', year{loop} ''',''fontsize'',12)'])
      ylabel('Y, km')
      xlabel('X, km')
      set(gca,'fontsize',12)
      eval(['print -dpng ' direc '/my4_locate3_workingfiltered' datestr(startDate(loop)+day-1,'yyyy-mm-dd') '.png'])
      %savefig([direc '/my_locate3_working' datestr(startDate(loop)+day-1,'yyyy-mm-dd')])
      %eval(['print -depsc ' direc '/my_locate3_working' datestr(startDate(loop)+day-1,'yyyy-mm-dd') '.eps'])
    end  
  end
  end
end