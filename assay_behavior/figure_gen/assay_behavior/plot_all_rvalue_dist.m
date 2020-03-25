% plot the r-value and p-value distributions for all collapsed representations
% of the decathlon data sets

% set decathlon_final_data.mat directory and load data
fdir = '';
D = load_decathlon_structs(fdir,'D13_als');

figure;
nreps = 20;
for i=1:2
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
    opts = {'CollapseFields';collapse_fields;'CollapseMode';collapse_mode;...
        'PCs';2;'Trim';true;'ImputeMode';'none'};
    D_p = D;
    data = cell(nreps,numel(D));
    if i==2
        for j=1:nreps
            D_p = D;
            D_p(1).data = D_p(1).data(randperm(size(D_p(1).data,1),192),:);
            D_p = pair_decathlon_structs(D_p,opts{:});
            for k=1:numel(D_p)
                d = D_p(k).data(~all(isnan(D_p(k).data),2),:);
                data{j,k} = d(randi(size(d,1),[192 1]),:);
            end
        end
    else
        D_p = pair_decathlon_structs(D_p,opts{:});
        for k=1:numel(D_p)
            d = D_p(k).data(~all(isnan(D_p(k).data),2),:);
            data(:,k) = arrayfun(@(ii) d(randi(size(d,1),[192 1]),:), 1:nreps,...
                'UniformOutput', false);
        end
    end

    % plot cross validation results
    ahs(1) = subplot(4,2,(i-1)*2+1);
    ahs(2) = subplot(4,2,(i-1)*2+2);
    [~,pvals,~,all_p]=plot_rvalue_dist(data,ahs);
    title(ahs(1),title_str);
    title(ahs(2),title_str);
    drawnow;
    
    bins = linspace(0,1,500);
    fdr_curves = cat(1,pvals{:,2})./cat(1,pvals{:,1});
    [~,alpha_idx] = min(abs(bins-0.05));
    for j=1:size(fdr_curves,1)
        fprintf('D%i FDR = %0.2f\n',j,fdr_curves(j,alpha_idx));
    end
end


% calculate the observed and expected overlapping significant p-vals
all_p = cat(2,all_p{:});
frac_sig = NaN(size(all_p,1),2);
for i=1:2   
   tmp_p = cellfun(@(p1,p2) cat(2,p1,p2), ...
       all_p(:,(i-1)*2+1), all_p(:,(i-1)*2+2),'UniformOutput', false);
   frac_sig(:,i) = cellfun(@(p) sum(all(p<0.05,2))/size(p,1), tmp_p);
end

