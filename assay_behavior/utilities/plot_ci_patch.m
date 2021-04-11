function [vx,vy] = plot_ci_patch(x, yvals, ci, color)

[vx,vy] = ci_vertices(x,yvals,ci); 
patch('XData',vx(:),'YData',vy(:),'FaceColor',color,'FaceAlpha',0.5,'EdgeColor','none');