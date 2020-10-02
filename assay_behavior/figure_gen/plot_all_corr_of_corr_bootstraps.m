% Generate r-value distribution for correlation of the decathlon correlation
% matrices to themselves and each other via bootstrapping raw data

% define data file path and load decathlon data
fdir = '';
D = load_decathlon_structs(fdir,'D13_als');

% iterate over different collapsed matrix representations
figure('Name','correlation of decathlon r-values');
n_modes = 2;
num_reps = 50;
for i=1:n_modes
    switch i
        case 1
            collapse_mode = 'average';
            collapse_fields = 'none';
            title_str = 'circ expanded';
        case 2
            collapse_mode = 'PCA';
            collapse_fields = 'all';
            title_str = 'a priori PCA';
    end
    opts = {'CollapseFields';collapse_fields;'CollapseMode';collapse_mode;'Trim';true};
    D_col = collapseMetrics(D,opts{:});
    
    % bootstrap data and plot distributions for observed data
    subplot(2,n_modes,i);
    bootstrap_dec_rval_corr(D_col,num_reps,mod(i,n_modes)>0,false);
    title(sprintf('observed r-value bootstrap (%s)',title_str));
    set(gca,'XLim',[-1 1],'XTick',-1:0.2:1);
    drawnow;
    
    % bootstrap data and plot distributions for shuffled data
    subplot(2,n_modes,i+n_modes);
    bootstrap_dec_rval_corr(D_col,num_reps,mod(i,n_modes)>0,true);
    title(sprintf('shuffled r-value bootstrap (%s)',title_str));
    set(gca,'XLim',[-1 1],'XTick',-1:0.2:1);
    drawnow;
end

