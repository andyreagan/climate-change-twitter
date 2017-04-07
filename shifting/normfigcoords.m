function [X,Y] = normfigcoords(x,y)
% [X,Y] = normfigcoords(x,y)
% 
% given normalized figure coordinates,
% returns in terms of current figures axes
%
% for use in general plots
%
% 0 < x < 1, 0 < y < 1
% 
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');

if (strcmp(get(gca,'xscale'),'linear'))
  X = x*(xlim(2)-xlim(1)) + xlim(1);
else
  X = xlim(2)^x/xlim(1)^(x-1);
end

if (strcmp(get(gca,'yscale'),'linear'))
  Y = y*(ylim(2)-ylim(1)) + ylim(1);
else
  Y = ylim(2)^y/ylim(1)^(y-1);
end

