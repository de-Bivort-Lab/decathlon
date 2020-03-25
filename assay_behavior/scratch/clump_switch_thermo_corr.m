
d = D_thermo.clumpiness.data;

[r,p] = cellfun(@(d) corr(d','rows','pairwise','Type','spearman'), d, 'UniformOutput', false);
r = cellfun(@(d) d(1,3), d);

figure;
histogram(r);

%%


c1 = cellfun(@(d,n) nansum(d(1,:).*n(1,:))/nansum(n(1,:)), d, n);
c2 = cellfun(@(d,n) nansum(d(idx2,:).*n(3,:))/nansum(n(3,:)), d, n);
f = ~isnan(c1) & ~isnan(c2);
c1 = c1(f);
c2 = c2(f);

r = corr(c1,c2);

%%

d = D_thermo.clumpiness.data;
s = D_thermo.switchiness.data;
idx1 = 3;
idx2 = 4;

figure;
ah=subplot(2,2,1);
c1 = nanzscore(cellfun(@(d,n) nanmean(d(idx1,:)), d));
c2 = nanzscore(cellfun(@(d,n) nanmean(d(idx2,:)), d));
f = ~isnan(c1) & ~isnan(c2);
c1 = c1(f);
c2 = c2(f);

r = corr(c1,c2,'Type','spearman');
PCARegressionCI([c1,c2],ah);
title(sprintf('r = %0.2f',r));
xlabel('clump (23C)');
ylabel('clump (29C)');

ah=subplot(2,2,2);
c1 = nanzscore(cellfun(@(d,n) nanmean(d(idx1,:)), s));
c2 = nanzscore(cellfun(@(d,n) nanmean(d(idx2,:)), s));
f = ~isnan(c1) & ~isnan(c2);
c1 = c1(f);
c2 = c2(f);

r = corr(c1,c2,'Type','spearman');
PCARegressionCI([c1,c2],ah);
title(sprintf('r = %0.2f',r));
xlabel('switch (23C)');
ylabel('switch (29C)');

ah=subplot(2,2,3);
c1 = nanzscore(cellfun(@(d,n) nanmean(d(idx1,:)), d));
c2 = nanzscore(cellfun(@(d,n) nanmean(d(idx1,:)), s));
f = ~isnan(c1) & ~isnan(c2);
c1 = c1(f);
c2 = c2(f);

r = corr(c1,c2,'Type','spearman');
PCARegressionCI([c1,c2],ah);
title(sprintf('r = %0.2f',r));
xlabel('clump (23C)');
ylabel('switch (23C)');

ah=subplot(2,2,4);
c1 = nanzscore(cellfun(@(d,n) nanmean(d(idx2,:)), d));
c2 = nanzscore(cellfun(@(d,n) nanmean(d(idx2,:)), s));
f = ~isnan(c1) & ~isnan(c2);
c1 = c1(f);
c2 = c2(f);

r = corr(c1,c2,'Type','spearman');
PCARegressionCI([c1,c2],ah);
title(sprintf('r = %0.2f',r));
xlabel('clump (29C)');
ylabel('switch (29C)');

%%

eff = D_thermo.effector;
effs = {'Shi';'Iso';'Trp'};
c = D_thermo.clumpiness.data;
s = D_thermo.switchiness.data;
colors = {[1 .5 0];[.5 .5 .5];[0 0 .8]};

figure;
uidx = unique_idx_pairs(3,1); 
ah = subplot_array(9);
for i=1:numel(ah)
    hold(ah(i),'on');
end

for j=1:size(uidx,1)
    
    lh = gobjects(numel(effs),3);
    r_labels = NaN(numel(effs),1);
    for i=1:numel(effs)
    
        curr_eff = strcmp(eff,effs{i});

        eff_c = cellfun(@(c) nanmean(c(1:3,:),2), c(curr_eff), 'UniformOutput', false);
        eff_c = nanzscore(cat(2,eff_c{:})');
        eff_s = cellfun(@(s) nanmean(s(1:3,:),2), s(curr_eff), 'UniformOutput', false);
        eff_s = nanzscore(cat(2,eff_s{:})');
        rc = corr(eff_c,'rows','pairwise','Type','spearman');
        rs = corr(eff_s,'rows','pairwise','Type','spearman');
        rcs = arrayfun(@(i) corr(eff_c(:,i),eff_s(:,i),'rows','pairwise','Type','spearman'),1:3);
    

        axes(ah((j-1)*3+1));
        r_labels(i,1) = rc(uidx(j,1),uidx(j,2));
        lh(i,1) = pretty_scatter(eff_c(:,uidx(j,1)),eff_c(:,uidx(j,2)),colors{i},'MarkerSize',3);
        xlabel(sprintf('clumpiness (block %i)',uidx(j,1)));
        ylabel(sprintf('clumpiness (block %i)',uidx(j,2)));
        
        axes(ah((j-1)*3+2));
        r_labels(i,2) = rs(uidx(j,1),uidx(j,2));
        lh(i,2) = pretty_scatter(eff_s(:,uidx(j,1)),eff_s(:,uidx(j,2)),colors{i},'MarkerSize',3);
        xlabel(sprintf('switchiness (block %i)',uidx(j,1)));
        ylabel(sprintf('switchiness (block %i)',uidx(j,2)));

        axes(ah((j-1)*3+3));
        r_labels(i,3) = rcs(j);
        lh(i,3) = pretty_scatter(eff_c(:,uidx(j,1)),eff_s(:,uidx(j,2)),colors{i},'MarkerSize',3);
        xlabel(sprintf('clumpiness (block %i)',j));
        ylabel(sprintf('switchiness (block %i)',j));
    end
    
    legend(lh(:,1),cellfun(@(e,r) sprintf('%s r=%0.2f',e,r),...
        effs, num2cell(r_labels(:,1)), 'UniformOutput', false));
    legend(lh(:,2),cellfun(@(e,r) sprintf('%s r=%0.2f',e,r),...
        effs, num2cell(r_labels(:,2)), 'UniformOutput', false));
    legend(lh(:,3),cellfun(@(e,r) sprintf('%s r=%0.2f',e,r),...
        effs, num2cell(r_labels(:,3)), 'UniformOutput', false));
end

