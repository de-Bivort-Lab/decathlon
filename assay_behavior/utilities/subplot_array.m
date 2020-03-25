function [ax_handles,n_rows,n_cols] = subplot_array(n_plots)

% calculate number of rows and columns in subplot grid (bias to more columns)
n_cols = ceil(sqrt(n_plots));
n_rows = ceil(n_plots/n_cols);

% initiatilze axes and store handles
ax_handles = gobjects(n_plots,1);
for i=1:numel(ax_handles)
    ax_handles(i) = subplot(n_rows,n_cols,i);
end
