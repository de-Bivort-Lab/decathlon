function [ax_handles,n_rows,n_cols] = subplot_array(n_plots,varargin)


ax_opts = {};
if ~isempty(varargin)
    n_rows = n_plots;
    n_cols = varargin{1};
    if numel(varargin) > 1
        ax_opts = varargin(2:end);
    end
else
    % calculate number of rows and columns in subplot grid (bias to more columns)
    n_cols = ceil(sqrt(n_plots));
    n_rows = ceil(n_plots/n_cols);
end

% initiatilze axes and store handles
ax_handles = gobjects(n_rows*n_cols,1);
for i=1:numel(ax_handles)
    ax_handles(i) = subplot(n_rows,n_cols,i,ax_opts{:});
end
