function line_handle = pretty_scatter(x,y,color,varargin)

% define default options
opts = {'Marker'; 'o'; 'LineStyle'; 'none'; 'MarkerEdgeColor'; 'none';...
        'MarkerSize'; 2; 'LineWidth'; 2};

% scatter data
line_handle = plot(x, y, 'MarkerFaceColor', color, opts{:},varargin{:});