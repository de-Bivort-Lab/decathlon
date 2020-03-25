function [corr_fh, pval_fh, r_vals] = plot_paired_corrmats(D,varargin)
% D = array of decathlon structs

% parse inputs
collapse_mode = 'average';
collapse_fields = 'circadian';
matches_only = false;
do_trim = true;
cluster = true;
pcs = 2;
for i=1:numel(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'CollapseFields'
                i = i+1;
                collapse_fields = varargin{i};
            case 'CollapseMode'
                i = i+1;
                collapse_mode = varargin{i};
            case 'MatchesOnly'
                i = i+1;
                matches_only = varargin{i};
            case 'PCs'
                i = i+1;
                pcs = varargin{i};
            case 'Trim'
                i = i+1;
                do_trim = varargin{i};
            case 'Cluster'
                i = i+1;
                cluster = varargin{i};
        end
    end
end
opts = parse_processing_options(varargin{:});

if strcmpi(collapse_fields,'all') || strcmpi(collapse_mode,'PCA')
   do_trim = false;
   trim_idx = find(strcmp(opts,'Trim'))+1;
   opts{trim_idx} = false;
end

% collapse data
if strcmpi(collapse_mode,'PCA')
    D = combine_PCA_split(D,opts{:});
else
    for i=1:numel(D)
        D(i) = collapseMetrics(D(i),opts{:});
    end
end

% initialize struct array for paired data
D_p = D;

% match fields between decathlon structs
for i=1:numel(D_p)
    for j=1:numel(D_p)-1
        other_idx = find(~ismember(1:numel(D_p),i));
        [D_p(i),D_p(other_idx(j))] = pairFields(D_p(i),D_p(other_idx(j)),opts{:});
    end
end

% find additional fields
D_u = D;
unique_fields = cell(numel(D),1);
if ~matches_only
    paired_fields = standardize_fieldnames(D_p(1).fields);
    if ~strcmpi(collapse_fields,'all')
        [assay, metric] = parse_fieldnames(paired_fields);
        paired_fields = cellfun(@(a,m) ...
            sprintf('%s %s',a,m), assay, metric, 'UniformOutput', false);
    end
    for i=1:numel(D)
        all_fields = standardize_fieldnames(D(i).fields);
        if ~strcmpi(collapse_fields,'all')
            [assay, metric, day] = parse_fieldnames(all_fields);
            all_fields = cellfun(@(a,m) ...
                sprintf('%s %s',a,m), assay, metric, 'UniformOutput', false);
            iscirc = strcmpi(assay,'Circadian');
            all_fields(iscirc) = cellfun(@(m,d) sprintf('Circadian %s (%i)',m,d), ...
                metric(iscirc), num2cell(day(iscirc)), 'UniformOutput', false);
        end
        unique_fields{i} = find(~any(str_list_contains(all_fields,paired_fields),2));
        D_u(i).fields = all_fields(unique_fields{i});
        D_u(i).data = D(i).data(:,unique_fields{i});
    end
end
max_dim = numel(D_p(1).fields) + max(cellfun(@numel,unique_fields));

% get labels and clustering permutation
if cluster
    [p, labels] = get_cluster_permutation(D_p);
    if numel(D_u)==3 && numel(D_u(2).fields)>1 && numel(D_u(3).fields)>1
        [D_u(2),D_u(3)] = pairFields(D_u(2),D_u(3),opts{:});
        [p_u, ~] = get_cluster_permutation(D_u(2:3));
        D_u(2:3) = apply_cluster_permutation(p_u, D_u(2:3));
    end
else
    p = 1:numel(D_p(1).fields);
    labels = D_p(1).fields;
end




% plot all matrices together
corr_fh = figure('Name','Sorted - paired correlation matrices');
np = numel(D_p(1).fields);
r_vals = cell(numel(D_p),1);
for i=1:numel(D_p)
    ah1 = subplot(2,numel(D_p),i);
    ah2 = subplot(2,numel(D_p),numel(D_p)+i);
    plot_dat = [D_p(i).data(:,p) D_u(i).data];
    plot_dat = [plot_dat NaN(size(plot_dat,1), max_dim-size(plot_dat,2))];
    plot_labels = [labels; D_u(i).fields];
    plot_labels = [plot_labels; repmat({''}, max_dim-numel(plot_labels), 1)];
    [~,r_vals{i}] = plotCorr(plot_dat, 'Labels', plot_labels, 'Cluster', false,...
       'Parent',[ah1; ah2],'Title',sprintf('decathlon %i',i), 'Patch', false,'Type','spearman');
    vx = [0 0 np np 0] + 0.5;
    vy = [0 np np 0 0] + 0.5;
    patch('XData',vx','YData',vy','FaceColor','none','EdgeColor','w',...
        'LineWidth',0.5, 'Parent', ah1);
end

all_ax = findall(corr_fh,'Type','axes');
for i=1:numel(all_ax)
    axis(all_ax(i),'equal','tight');
end

rh = groot;
pos = rh.MonitorPositions(end,:).*[1.01 0 0.8 0.8];
pos(2) = pos(2) - 50;
corr_fh.Units = 'pixels';
corr_fh.Position = pos;
pval_fh.Units = 'pixels';
pval_fh.Position = pos;