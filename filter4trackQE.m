function [locTrue, locFalse] = filter4track(loc,drMax,dtMax,n4track)
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



x = [loc.x];
y = [loc.y];
ot = [loc.ot];
n = length(x);
keep = false(n,1);
throw = false(n,1);


for i=1:n
  dr = sqrt((x-x(i)).^2+(y-y(i)).^2);
  dt = abs(ot - ot(i));
  near = dr<drMax & dt<dtMax;
  if sum(near)>n4track
    keep(near) = true;

  end
end
[indices]=find(keep == 0);


locTrue = loc(keep);
locFalse = loc(indices);