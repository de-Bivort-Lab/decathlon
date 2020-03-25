function plot_pair_dot(a,b,color,patch_opts)

% define default options
opts = {'Marker'; 'o'; 'LineStyle'; 'none'; 'MarkerEdgeColor'; 'none';...
        'MarkerSize'; 2; 'LineWidth'; 2};
if isempty(patch_opts)
    patch_opts = {};
end
    

% initialize position vectors
x = normrnd(1,0.05,[numel(a) 2]);
x(:,2) = x(:,2)+1;
y = [a b];

% plot points
plot(x(:),y(:),'MarkerFaceColor',color,opts{:});
min_y = min(y(:));
max_y = max(y(:));
range_y = max_y-min_y;
set(gca,'XLim',[0.5 2.5],'YLim',[min_y-(range_y)*0.1 max_y+(range_y)*0.1]);

patch('Faces',1:numel(a),'XData',x','YData',y','FaceColor','none',...
    'EdgeColor',color,'EdgeAlpha',0.1,patch_opts{:});