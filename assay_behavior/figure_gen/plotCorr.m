function [varargout] = plotCorr(data,varargin)
% Plot pairwise correlation matrices and (optionally) sort rows and columns
% via hierarchical clustering. Additionally plot accompanying
% pairwise correlation p-value matrices

% parse inputs
fs = 6;
cluster = true;
save_path = '';
plot_title = '';
signed = true;
ext = {'.fig'};
alpha = 0.01;
corr_type = 'pearson';
corr_opts = {'rows';'pairwise'};
pval_patch = true;
parent = [];
plot_pvalues = false;
for i=1:length(varargin)
    
    arg = varargin{i};
    if ischar(arg)
    switch arg
        case 'Labels'
            i=i+1;
            labels = varargin{i};       % column labels for raw data
        case 'FontSize'
            i=i+1;
            fs = varargin{i};           % plot font size
        case 'Cluster'
            i=i+1;
            cluster = varargin{i};
        case 'SavePath'
            i = i+1;
            save_path = varargin{i};
        case 'Title'
            i=i+1;
            plot_title = varargin{i};
        case 'Signed'
            i=i+1;
            signed = varargin{i};
        case 'Ext'
            i=i+1;
            ext = varargin{i};
        case 'Parent'
            i=i+1;
            parent = varargin{i};
        case 'PvalPatch'
            i=i+1;
            pval_patch = varargin{i};
        case 'Type'
            i = i+1;
            corr_type = varargin{i};
        case 'Options'
            i=i+1;
            corr_opts = varargin{i};
        case 'PvalPlot'
            i=i+1;
            plot_pvalues = varargin{i};
    end
    end
end

% calculate covariance matrix
[r,p] = corrcoef(data,corr_opts{:});

% replace NaNs
r(isnan(r))=0;
p(isnan(p))=1;

% sort rows and columns by hierarchical clustering
Zoutperm = 1:length(r);
if cluster
    data(all(isnan(data),2),:) = [];
    Z=linkage(data','average','correlation');
%     Z=linkage(r);
    f=figure;
    [ZH, ZT, Zoutperm]=dendrogram(Z,0);
    close(f);
    r=r(Zoutperm,Zoutperm);
    p=p(Zoutperm,Zoutperm);
end

% plot correlation matrix
if ~isempty(parent)
    ah1 = parent(1); 
    fh1 = ah1.Parent;
else
    fh1 = figure;
    ah1 = axes;
end
if signed
    imh = imagesc(r,'Parent',ah1);
    caxis(ah1,[-1,1]);
else
    imh = imagesc(abs(r),'Parent',ah1);
    caxis(ah1,[0,1]);
end

% configure corr mat axes
colormap(ah1,nanticoke);
colorbar(ah1);
set(ah1,'TickLength',[0 0]);    

if ~isempty(plot_title)
   title(sprintf('%s (r-vals)',plot_title)); 
end

if exist('labels','var')
    clusteredLabels=labels(Zoutperm);

    % format field labels for display
    clusteredLabels = pretty_labels(clusteredLabels);
    
    set(ah1,'Ytick',1:length(clusteredLabels),'YtickLabel', clusteredLabels,'fontsize',fs);
    set(ah1,'XTick',1:length(labels),'XTickLabel',clusteredLabels,'fontsize',fs,'XTickLabelRotation',90);
else
    clusteredLabels = repmat({''},size(r,1),1);
end

if ~isempty(save_path)
    for i=1:numel(ext)
        saveas(fh1,sprintf('%s%s%s',save_path,'_corrmat',ext{i}));
    end
end

% plot pvalues
fh2 = gcf;
if plot_pvalues
    if pval_patch
       [row_idx,col_idx] = find(p < alpha);
       vx = repmat(col_idx',5,1) + repmat(0.5.*[-1 -1 1 1 -1]',1,numel(col_idx));
       vy = repmat(row_idx',5,1) + repmat(0.5.*[1 -1 -1 1 1]',1,numel(col_idx));
       patch('XData',vx,'YData',vy,'FaceColor','none','EdgeColor','w',...
           'LineWidth',0.25,'Parent',ah1);
    end

    if numel(parent)>1
        ah2 = parent(2); 
        fh2 = ah1.Parent;
    else
        fh2 = figure;
        ah2 = axes;
    end
    imagesc(p,'Parent',ah2);
    colorbar(ah2);
    c=[0 1 1];
    logcmap =interp1(fliplr(logspace(0,log10(256),5)),...
        [0 0 0; .25 .25 .25; 0 1 1; 0 0 1; 1 0 1], 1:256);
    colormap(ah2,logcmap);
    set(gca,'Ytick',1:length(clusteredLabels),'YtickLabel', clusteredLabels,'fontsize',fs);
    set(gca,'XTick',1:length(clusteredLabels),'XTickLabel',clusteredLabels,'fontsize',fs,'XTickLabelRotation',90);
    set(gca,'TickLength',[0 0]);

    if ~isempty(plot_title)
       title(sprintf('%s (p-vals)',plot_title)); 
    end

    if ~isempty(save_path)
        for i=1:numel(ext)
            saveas(fh2,sprintf('%s%s%s',save_path,'_pvals',ext{i}));
        end
    end
end

% parse outputs
for i = 1:nargout
    switch i
        case 1, varargout(i)={[fh1;fh2]};
        case 2, varargout(i)={r};
        case 3, varargout(i)={p};
        case 4, varargout(i)={Zoutperm};
    end
end

