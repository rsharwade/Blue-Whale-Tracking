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

minQEtrue=.15
Month='April';
year1='2012';
datefirst='04/01/2012';
datelast='04/30/2012';
date=[datenum(datefirst):1:datenum(datelast)];
date1=date(1)-1;
QEtrues=[];
QEfalses=[];
QEtruesSTD=[];
QEfalsesSTD=[];
day=[];
clf
for k=1:length(date)




loc2date=[];
for j=1:length(loc2);
    if floor(loc2(j).ot) == date(k);
        
        loc2date=[loc2date loc2(j)];
     
    end
end



  
if length(loc2date) > 0
  
  % Filter loc2 locations for tracks
  [locTrue, locFalse] = filter4trackQE(loc2date,drMax,dtMax,n4track);
  ot = [loc2.ot]-datenum(date(k))+1;
  x = [loc2.x];
  y = [loc2.y];
  otFilt = [locTrue.ot]-datenum(date(k))+1;
  xFilt = [locTrue.x];
  yFilt = [locTrue.y];
  xfalse = [locFalse.x];
  yfalse = [locFalse.y];
  xall = [loc2date.x];
  yall = [loc2date.y];

 %plot(staVec.x,staVec.y,'k^','MarkerSize',7,'LineWidth',2);
 %title([datestr(date(k)) 'filtered'])
%hold on
%plot(xall,yall,'mo')
 %plot(xFilt,yFilt,'ko')
% plot(xfalse,yfalse,'go')
 
 %ginput(1)
%clf

 QEtrue=[];
 for j=1:length(locTrue);
 
        QEtrue=[QEtrue locTrue(j).QEmax];
     
 end

  
 QEfalse=[];
 for j=1:length(locFalse);
 
        QEfalse=[QEfalse locFalse(j).QEmax];
     
    end
    

 
day1=date(k)-date1;  
meanQEfalse=mean(QEfalse);
meanQEtrue=mean(QEtrue);
stdQEfalse=std(QEfalse);
stdQEtrue=std(QEtrue);

QEtrues=[QEtrues meanQEtrue];
QEfalses=[QEfalses meanQEfalse];
QEtruesSTD=[QEtruesSTD stdQEtrue];
QEfalsesSTD=[QEfalsesSTD stdQEfalse];
day=[day day1];

else disp(['no data' datestr(date(k))])
end

end

QEvalues=[day.' QEtrues.' QEtruesSTD.' QEfalses.' QEfalsesSTD.'] ;
QEdiff=[QEvalues(:,2)-QEvalues(:,4)];

QEstrong=[];
QEweak=[];
for v=1:length(QEdiff(:,1));
    if QEdiff(v) > minQEtrue;
        QEstrong=[QEstrong; QEvalues(v,:) QEdiff(v)];
    else QEweak=[QEweak; QEvalues(v,:) QEdiff(v)];
end
end

QEvalues
QEstrong
QEweak


bar(QEvalues(:,1), QEvalues(:,2),'g')
hold on
bar(QEweak(:,1),QEweak(:,2),'y')
bar(QEvalues(:,1), QEvalues(:,4),'r')

ylabel('Average QE value')
xlabel('Days')
legend('Strong correct locations','Weak correct locations', 'False locations')
title([Month year1 ' QE comparison'])
saveas(gcf,[Month year1 'QEbar.jpg'],'jpg')