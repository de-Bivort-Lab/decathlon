D = load_decathlon_structs;
D = avg_als_impute(D,200);

%%

opts = {'CollapseMode','PCA','CollapseFields','all'};
D_p = pair_decathlon_structs(D,opts{:});

%%

save_dir = 'C:\Users\winsl0w\Documents\decathlon\decathlon_analysis\figures\pca_loadings\';

plot_pca_loadings(D_p(1),'SaveDir',save_dir);

fig_paths = recursiveSearch(save_dir,'ext','.fig');
for i=1:numel(fig_paths)
    fh = open(fig_paths{i});
    [dd,nn,~] = fileparts(fig_paths{i});
    saveas(fh,sprintf('%s/%s.pdf',dd,nn));
    close(fh);
end