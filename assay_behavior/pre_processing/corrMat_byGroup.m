function corrMat_byGroup(groups, D, varargin)

plot_corr_options = {};
for i=1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
    switch arg
        case 'PlotCorrOptions'
            i=i+1;
            plot_corr_options = varargin{i};
    end
    end
end

for i=1:numel(groups)
    
    [f, fidx] = groupFields(D.fields,groups{i});
    grp_mat = D.data(:,fidx);
    grp_label = groups{i};
    if iscell(grp_label)
        grp_label = grp_label{1};
    end
    plotCorr(grp_mat,'Labels',f,'Title',grp_label,plot_corr_options{:});
end