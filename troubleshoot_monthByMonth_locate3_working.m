% Script to run locations 1 month at a time using locate3_working and create daily plots
% Some plot limits are hardwired
% *** Note that for this to work 3 parameters need to be commented out in locate3_working ***

%% Parameters
% Criteria to group calls in a track
drMax = 10;   %Maximum allowed spacing in kilometers
dtMax = 1/24; %Maximum spacing in time (days)
n4track = 4;  %Number of other calls nearby in space and time to consider a call to be in a track

% Number of calls in a track for it to be signitifant (which means it gets plotted)
nSignificant = 10;

% directory for plots
direc = '/burr1/CI/locations/monthly_locations3';

% Files in wich to save solutions for each interval
file = {'monthly_locations3/new100location_Feb2012.mat'};
      
% Time intervals to process
month = {'February' 'March'};
year = {'2012' '2012'};
startDate = [datenum('Feb 1 2012')];
endDate = [datenum('Mar 1 2012')];

%% Loop over intervals to process
for loop = 1:length(file)
  
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

  %% Make daily plots
  for day = 6
    
    % Plot loc2 locations
    j = find(ot>=day & ot<day+1);
    if ~isempty(j) 
      figure(201); clf
      plot([staVec.x(1:end)],[staVec.y(1:end)],'sk','markerfacecolor','k','markersize',9)
      for i=[40 42 36 31 23 4 12 16 13]
        text(staVec.x(i),staVec.y(i),staVec.name(i),'fontsize',12)
      end
      hold on
      if length(j)>1
        col = LinColor((ot(j)-floor(min(ot(j(1)))))*24,[0 24],'gray');
      else
        col = [0 0 1];
      end
      for i=1:length(j)
        plot3(x(j(i)),y(j(i)),repmat(100,length(x(j(i))),1),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:));
      end
      axis equal
   
      % Hardwire plot limits so each day has the same scale
      xlim([min(staVec.x)-10 max(staVec.x)+50])
      ylim([min(staVec.y)-10 max(staVec.y)+10])
      eval(['title('' Locate3 Working - ' month{loop} ' ' int2str(day) ', ', year{loop} ''',''fontsize'',12)'])
      ylabel('Y, km')
      xlabel('X, km')
      eval(['print -dpng ' direc '/newlocate3_working_100km' datestr(startDate(loop)+day-1,'yyyy-mm-dd') '.png'])
      savefig([direc '/newlocate3_working_100km' datestr(startDate(loop)+day-1,'yyyy-mm-dd')])

%       eval(['print -depsc ' direc '\locate3_working_' month{loop} int2str(day) '.eps'])
    end
    
    % plot loc2 locations after filtering for tracks
    j = find(otFilt>=day & otFilt<day+1);
    if ~isempty(j) && length(j)>nSignificant
      figure(1); %clf
      plot([staVec.x(1:end)],[staVec.y(1:end)],'sk','markerfacecolor','k','markersize',9)
      for i=[40 42 36 31 23 4 12 16 13]
        text(staVec.x(i),staVec.y(i),staVec.name(i),'fontsize',12)
      end
      hold on
      if length(j)>1
        col = LinColor((otFilt(j)-floor(min(otFilt(j(1)))))*24,[0 24],'gray');
      else
        col = [0 0 1];
      end
      for i=1:length(j)
        plot3(xFilt(j(i)),yFilt(j(i)),repmat(100,length(xFilt(j(i))),1),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:));
      end
      axis equal
    
      % Hardwire plot limits so each day has the same scale
      xlim([min(staVec.x)-10 max(staVec.x)+50])
      ylim([min(staVec.y)-10 max(staVec.y)+10])
      eval(['title('' Locate3 Working Filtered - ' month{loop} ' ' int2str(day) ', ', year{loop} ''',''fontsize'',12)'])
      ylabel('Y, km')
      xlabel('X, km')
      set(gca,'fontsize',12)
      eval(['print -dpng ' direc '/newfiltered_locate3_working_100km' datestr(startDate(loop)+day-1,'yyyy-mm-dd') '.png'])
      savefig([direc '/newfiltered_locate3_working_100km' datestr(startDate(loop)+day-1,'yyyy-mm-dd')])
      eval(['print -depsc ' direc '\newfiltered_locate3_working_100km' datestr(startDate(loop)+day-1,'yyyy-mm-dd') '.eps'])
    end  
  end
  end
end