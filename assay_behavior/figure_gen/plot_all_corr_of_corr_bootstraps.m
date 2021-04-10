% Generate r-value distribution for correlation of the decathlon correlation
% matrices to themselves and each other via bootstrapping raw data

% define data file path and load decathlon data
D = load_decathlon_structs(pwd,'D_als_filled');

%% iterate over different collapsed matrix representations
figure('Name','correlation of decathlon r-values');
n_modes = 2;
num_reps = 50;
for i=1:n_modes
    switch i
        case 1
            collapse_mode = 'average';
            collapse_fields = 'none';
            title_str = 'Full';
        case 2
            collapse_mode = 'PCA';
            collapse_fields = 'all';
            title_str = 'Distilled';
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

%% Plot R-value scatter plots

pairs = unique_idx_pairs(numel(D),1);
nplots = size(pairs,1) * n_modes;
f2 = figure;
for i=1:n_modes
    switch i
        case 1
            collapse_mode = 'average';
            collapse_fields = 'none';
            title_str = 'Full';
        case 2
            collapse_mode = 'PCA';
            collapse_fields = 'all';
            title_str = 'Distilled';
    end
    opts = {'CollapseFields';collapse_fields;'CollapseMode';collapse_mode;'Trim';true};
    D_p = pair_decathlon_structs(D,opts{:});
    
    
    for j=1:size(pairs,1)
        ah = subplot(n_modes,size(pairs,1),(i-1)*size(pairs,1) + j);
       [r,~] = corr_of_corrcoef(D_p(pairs(j,1)), D_p(pairs(j,2)), 'Plot', true, 'Title', title_str);
       xlabel(sprintf('decathlon %i',pairs(j,1)));
       ylabel(sprintf('decathlon %i',pairs(j,2)));
    end
end

ahs = findall(f2,'Type','axes');
set(ahs,'Units','inches');
for i=1:numel(ahs)
   ahs(i).Position(3:4) = 0.75; 
end



