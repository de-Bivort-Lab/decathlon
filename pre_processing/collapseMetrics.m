function D = collapseMetrics(D,varargin)



% parse inputs
mode = 'average';
fields = 'circadian';
PCs = 2;
grp_path = 'D:\decathlon_data_and_analysis\decathlon 8-2017\meta\apriori_groups.mat';
for i=1:length(varargin)
    
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'CollapseMode'
                i=i+1;
                mode = varargin{i};
            case 'CollapseFields'
                i=i+1;
                fields = varargin{i};
            case 'PCs'
                i=i+1;
                PCs = varargin{i};
        end
    end
end

% iterate recursively if more than one element
if numel(D)>1
    if strcmpi(mode,'PCA')
        D = combine_PCA_split(D,varargin{:});
    else
        for i=1:numel(D)
            D(i) = collapseMetrics(D(i),varargin{:});
        end
    end
    return
end

f = D.fields;
    


switch lower(fields)
    case 'circadian'
        
        % find relevant indices in the matrix
        idx = find(~cellfun(@isempty,strfind(f,'Circadian')));
        cf = f(idx);
        cf = cellfun(@(x,y) x(1:y-2),cf,strfind(cf,'('),'UniformOutput',false);
        cf = cellfun(@(x,y) x(y+10:end),cf,strfind(cf,'Circadian'),'UniformOutput',false);
        
        % get unique field names and compute average
        uf = unique(cf);
        ufilt = cellfun(@(x) strcmp(cf,x),uf,'UniformOutput',false);
        collapsed = cellfun(@(x) nanmean(D.data(:,idx(x)),2),ufilt,'UniformOutput',false);
        newFields = cellfun(@(x) ['Circadian ' x],uf,'UniformOutput',false);
        
        % remove raw data and old field names
        D.data(:,idx)=[];
        D.fields(idx)=[];
        D.data = [D.data cat(2,collapsed{:})];
        D.fields = [D.fields; newFields];      
        
    case 'all'
        
        [apriori_data, apriori_names, grp_idx] = group_apriori_fields(D);
        empty_grp = cellfun(@isempty,grp_idx);
        apriori_data(empty_grp) = [];
        apriori_names(empty_grp) = [];
        grp_idx(empty_grp) = [];

        switch mode
            case 'average'
                collapsed_data = cellfun(@(ad) nanmean(ad,2), ...
                            apriori_data, 'UniformOutput', false);
                dMat = cat(2,collapsed_data{:});
                nf = apriori_names;
            case 'PCA'
                [~, ~, bs_exp] = cellfun(@(d) bootstrap_pca_nullmodel(d,100,size(d,1)),...
                    apriori_data,'UniformOutput',false);
                bs_exp = cellfun(@(b) cat(2,b{:}), bs_exp, 'UniformOutput', false);
                low_ci_bound = cellfun(@(b) prctile(b,2.5,2), bs_exp, 'UniformOutput', false);
                [coef,score,lat,~,explained] = ...
                    cellfun(@pca,apriori_data,'UniformOutput',false);
   
                % determine number of PCs to keep
                npcs = cellfun(@(ex,ci) find(ex>=ci,1,'Last'), explained, low_ci_bound,'UniformOutput', false);
                npcs(cellfun(@isempty,npcs)) = {1};
                npcs = cat(1,npcs{:});
                
                % plot PCA bootstrap
%                 figure;
%                 ncol = ceil(sqrt(numel(npcs)));
%                 nrow = ceil(numel(npcs)/ncol);
%                 for i=1:numel(npcs)
%                     subplot(nrow,ncol,i);
%                     plot_pca_bootstrap(apriori_data{i},100,95,'noncummulative',[1 0 0]);
%                     set(gca,'XLim',[0.99 size(apriori_data{i},2)],'YLim',[0 100]);
%                     title(sprintf('%s, k=%i',apriori_names{i},npcs(i)));
%                 end
                
                pc_data = cellfun(@(sc,n) sc(:,1:n), score, ...
                            num2cell(npcs),'UniformOutput', false);
                %dMat = nanzscore(cat(2,pc_data{:}));
                dMat = cat(2,pc_data{:});
                % create new field names
                nf = cell(sum(npcs),1);
                j = 0;
                for i=1:numel(apriori_names)
                    nf(j+1:j+npcs(i)) = arrayfun(@(ii) ...
                        sprintf('%s (PC%i)',apriori_names{i},ii),...
                        1:npcs(i), 'UniformOutput', false);
                    j = j+npcs(i);
                end

                D.loadings = cellfun(@(c,n) num2cell(c(:,1:n),1), ...
                    coef, num2cell(npcs), 'UniformOutput', false);
                D.loadings = cat(2,D.loadings{:})';
                D.loadings_labels = cellfun(@(gi,n) repmat({D.fields(gi)},n,1),...
                    grp_idx, num2cell(npcs), 'UniformOutput', false);
                D.loadings_labels = cat(1,D.loadings_labels{:});
        end
        
        D.data = dMat;
        D.fields = nf;
end

