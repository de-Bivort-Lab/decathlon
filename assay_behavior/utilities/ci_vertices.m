function [vx,vy] = ci_vertices(x, yvals, ci)

lower = prctile(yvals,(100-ci)/2);
upper = prctile(yvals,100-(100-ci)/2);

vx = [x(1) x fliplr(x)];
vy = [upper(1) lower  fliplr(upper)];