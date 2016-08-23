function loc = filter4track(loc,drMax,dtMax,n4track)
% Filter a structure of locations for tracks based on groups of call that are nearby in space and time
%
% Usage
%   loc = filter4track(loc,drMax,dtMax,n4track);
%
% Inputs
%   loc     - Location structure
%   drMax   - Maximum spacing in distance (km) to group calls
%   dtMax   - Maximum spacing in time (days) to group calls
%   n4track - If a call is near this many other calls in space and time then it is in a track
%
% Outputs
%   loc - Location structure that is limited to locations that are in a track

calls=[];
locnumber=[];

rmstry=[];
for j=1:length(loc); 
    rmstry=[rmstry loc(j).rms]; 
end

[rmslist,rmsinds]=sort(rmstry,'descend');

loc=loc(rmsinds);


for j=1:length(loc);
    calls=[calls; loc(j).det.peakAmp(loc(j).det.inLocation).*loc(j).det.duration(loc(j).det.inLocation)]; 
end

for j=1:length(loc);
    locnumber=[locnumber; j.*(ones(length(loc(j).det.duration(loc(j).det.inLocation)),1))];
end

[strengths,inds]=unique(calls);
a=1
usedlocs=[];
for q=1:length(strengths)
    
    locinds=find(calls == strengths(q));
    if length(locinds) == 1
        if isempty(find(usedlocs == locnumber(locinds))) == 1;
        loc1(a)=loc(locnumber(locinds));
        usedlocs=[usedlocs locnumber(locinds)];
        a=a+1;
        end
    end
    
    [bestval,minind]=min([loc(locnumber(locinds)).rms]);
    if isempty(find(usedlocs == locnumber(locinds(minind)))) == 1;
    loc1(a)=loc(locnumber(locinds(minind)));
    usedlocs=[usedlocs locnumber(locinds(minind))];
    a=a+1;
    end
end
timeloc=[];

for j=1:length(loc1); 
    timeloc=[timeloc loc1(j).ot]; 
end

[timelist,timeinds]=sort(timeloc,'ascend');

loc=loc1(timeinds);



x = [loc.x];
y = [loc.y];
ot = [loc.ot];
n = length(x);
keep = false(n,1);

for i=1:n
    dr = sqrt((x-x(i)).^2+(y-y(i)).^2);
    dt = abs(ot - ot(i));
  near = dr<drMax & dt<dtMax;
  if sum(near)>n4track
    keep(near) = true;
  end
end

loc = loc(keep);