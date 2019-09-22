% Generates plots from the descending neuron (DN) screen behavioral 
% data (Cande et. al, eLife 2018). 
% Load the DN data matrices (DN_ctl and DN_exp) before running.


% Plot correlation matrices and scree plots
fh1 = figure;
fh2 = figure;
effective_dims = NaN(4,1);
labels = {'Exp. light OFF';'Exp. light ON';'Ctl. light OFF';'Ctl. light ON'};
for i=1:4
    switch i
        case 1 
            %data = DN_exp.line_avg_off_pdfs;
            data = DN_exp.dn_avg_off_pdfs;
            cluster = true;
            p=1:size(DN_ctl.line_avg_off_pdfs,2);
        case 2 
            %data = DN_exp.line_avg_on_pdfs;
            data = DN_exp.dn_avg_on_pdfs;
            cluster = false;
            p=p_out;
        case 3 
            %data = DN_ctl.line_avg_off_pdfs;
            data = DN_ctl.dn_avg_off_pdfs;
            cluster = false;
        case 4
            %data = DN_ctl.line_avg_on_pdfs;
            data = DN_ctl.dn_avg_on_pdfs;
            cluster = false;
    end
    
    figure(fh1);
    ah = subplot(2,2,i);
    [~,~,~,p_out]=plotCorr(data(:,p),'Patch',false,'Parent',ah,'Cluster',cluster);
    axis(ah,'equal','tight');
    title(ah,labels{i});
    
    figure(fh2);
    ah = subplot(2,2,i);
    [~,~,nkeep] = plot_pca_bootstrap(data,100,95,'noncummulative',[.7 .7 .7]);
    set(gca,'YLim',[1E-2 1E2]);
    title(ah,{labels{i};sprintf('thresh = %i',nkeep)});
    effective_dims(i) = nkeep;
end

%% Compute the effective dimensionality of individual behaviors in each line separately

k_effective = NaN(numel(unique_DNs),4);
fh = figure;
for i=1:numel(unique_DNs)
    
    fprintf('iteration %i of %i\n',i,numel(unique_DNs));
    
    d = exp_off_pdfs{i}(~any(isnan(exp_off_pdfs{i}),2),:);
    [~,~,k_effective(i,1),ph] = plot_pca_bootstrap(d,100,95,'noncummulative',[.7 .7 .7]);
    cellfun(@delete,ph);
    
    d = exp_on_pdfs{i}(~any(isnan(exp_on_pdfs{i}),2),:);
    [~,~,k_effective(i,2),ph] = plot_pca_bootstrap(d,100,95,'noncummulative',[.7 .7 .7]);
    cellfun(@delete,ph);

    d = ctl_off_pdfs{i}(~any(isnan(ctl_off_pdfs{i}),2),:);
    [~,~,k_effective(i,3),ph] = plot_pca_bootstrap(d,100,95,'noncummulative',[.7 .7 .7]);
    cellfun(@delete,ph);

    d = ctl_on_pdfs{i}(~any(isnan(ctl_on_pdfs{i}),2),:);
    [~,~,k_effective(i,4),ph] = plot_pca_bootstrap(d,100,95,'noncummulative',[.7 .7 .7]);
    cellfun(@delete,ph);
end
delete(fh);

figure;
k_effective = k_effective(~any(isnan(k_effective),2),:);
ci = [bootstrap_mean_CI(k_effective,0.05,100)';NaN(1,size(k_effective,2))];
vx = [repmat(1:size(k_effective,2),2,1);NaN(1,size(k_effective,2))];
bar(mean(k_effective)); hold on;
plot(vx(:),ci(:),'k-','LineWidth',1);
ylabel('effective dimensionality');
set(gca,'XTick',1:size(k_effective,2),'TickLength',[0 0],...
    'XtickLabel',labels,'Xticklabelrotation',90,'YLim',[0 5]);
title('DN screen individual dimensionality');



%% Descending Neuron Screen t-SNE plots

% concatenate data matrices
all_exp_data = zscore([DN_exp.dn_avg_off_pdfs;DN_exp.dn_avg_on_pdfs]);
exp_off_z_data = all_exp_data(1:size(DN_ctl.dn_avg_off_pdfs,1),:);
exp_on_z_data = all_exp_data(size(DN_ctl.dn_avg_on_pdfs,1)+1:end,:);

all_ctl_data = zscore([DN_ctl.dn_avg_off_pdfs;DN_ctl.dn_avg_on_pdfs]);
ctl_off_z_data = all_ctl_data(1:size(DN_ctl.dn_avg_off_pdfs,1),:);
ctl_on_z_data = all_ctl_data(size(DN_ctl.dn_avg_on_pdfs,1)+1:end,:);

% initialize settings
fh = figure;
dist_type = 'euclidean';
perplex = 30;
fdir = 'C:\Users\winsl0w\Documents\decathlon\decathlon_analysis\matrices\decathlon_paper\cande_data';
fname = sprintf('Cande_tsne_perplex%i_%s.fig',perplex,dist_type);
ncol = 2;
opts = statset;
opts.MaxIter = 5000;
opts.TolFun = 1E-12;
opts.Display = 'iter';

% tSNE individuals (experimental)
ah = subplot(2,ncol,1); hold on;
embedded_data_a = tsne(all_exp_data,'Distance',dist_type,'Perplexity',perplex,'Options',opts);
mask_a = false(size(all_exp_data,1),1);
mask_a(1:size(exp_on_z_data,1)) = true;
pretty_scatter(embedded_data_a(mask_a,1),embedded_data_a(mask_a,2),'k');
pretty_scatter(embedded_data_a(~mask_a,1),embedded_data_a(~mask_a,2),'r');
title({'Cande data (Experimentals)';sprintf('t-SNE %s - DNs',dist_type)});
axis('equal');
xlim = [min(embedded_data_a(:,1)) max(embedded_data_a(:,1))];
ylim = [min(embedded_data_a(:,2)) max(embedded_data_a(:,2))];
buff = max([diff(xlim) diff(ylim)])*0.57;
xlim = [mean(xlim)-buff mean(xlim)+buff];
ylim = [mean(ylim)-buff mean(ylim)+buff];
set(gca,'XLim',xlim,'YLim',ylim);
legend({'lights OFF','lights ON'});

% tSNE metrics (experimental)
ah = subplot(2,ncol,2); hold on;
embedded_data_b = tsne(all_exp_data','Distance',dist_type,...
    'Perplexity',perplex,'Algorithm','exact','Options',opts);
ah1=pretty_scatter(embedded_data_b(:,1),embedded_data_b(:,2),'k');
title({'Cande data (Experimentals)';sprintf('t-SNE %s - metrics',dist_type)});
axis('equal');
set(gca,'XLim',[-25 25],'YLim',[-25 25]);


ah = subplot(2,ncol,3); hold on;
embedded_data_a = tsne(all_ctl_data,'Distance',dist_type,'Perplexity',perplex,'Options',opts);
mask_a = false(size(all_ctl_data,1),1);
mask_a(1:size(ctl_on_z_data,1)) = true;
pretty_scatter(embedded_data_a(mask_a,1),embedded_data_a(mask_a,2),'k');
pretty_scatter(embedded_data_a(~mask_a,1),embedded_data_a(~mask_a,2),'r');
title({'Cande data (Controls)';sprintf('t-SNE %s - DNs',dist_type)});
axis('equal');
xlim = [min(embedded_data_a(:,1)) max(embedded_data_a(:,1))];
ylim = [min(embedded_data_a(:,2)) max(embedded_data_a(:,2))];
buff = max([diff(xlim) diff(ylim)])*0.57;
xlim = [mean(xlim)-buff mean(xlim)+buff];
ylim = [mean(ylim)-buff mean(ylim)+buff];
set(gca,'XLim',xlim,'YLim',ylim);
legend({'lights OFF','lights ON'});


% concatenate grouped data
ah = subplot(2,ncol,4); hold on;
embedded_data_b = tsne(all_ctl_data','Distance',dist_type,...
    'Perplexity',perplex,'Algorithm','exact','Options',opts);
ah1=pretty_scatter(embedded_data_b(:,1),embedded_data_b(:,2),'k');
title({'Cande data (Experimentals)';sprintf('t-SNE %s - metrics',dist_type)});
axis('equal');
set(gca,'XLim',[-25 25],'YLim',[-25 25]);

set(findall(fh,'Type','axes'),'XTick',[],'YTick',[],'TickLength',[0 0]);
savefig(fh,[fdir fname]);

%% Compute effective dimensionality (connected components of correlation matrix)

D = struct('data',[],'fields',{});
D(1).data = exp_off_z_data;
D(2).data = exp_on_z_data;
D(3).data = ctl_off_z_data;
D(4).data = ctl_on_z_data;
d_label = {'lights off';'lights on'};

figure; 
subplot(1,2,1); hold on;
lhs = gobjects(2,1);
colors = {[0 0 0];[1 0 0]};
for i=1:2
    d = D(i).data;
    [out,intrinsic_dim,~]=decathlonConnCompSweep(d,1000);
    bins = 1:2:size(d,2);
    cts = histc(intrinsic_dim,bins);
    lhs(i) = bar(bins,log(cts),'EdgeColor','none',...
        'FaceColor',colors{i},'FaceAlpha',0.5);
    ylabel('log(count)');
    xlabel('effective dimensionality');
end
set(gca,'XLim',[0 size(d,2)]);
legend(lhs,d_label);
title('DN screen (experimentals)');

subplot(1,2,2); hold on;
lhs = gobjects(2,1);
colors = {[0 0 0];[1 0 0]};
for i=1:2
    d = D(i+2).data;
    [out,intrinsic_dim,~]=decathlonConnCompSweep(d,1000);
    bins = 1:2:size(d,2);
    cts = histc(intrinsic_dim,bins);
    lhs(i) = bar(bins,log(cts),'EdgeColor','none',...
        'FaceColor',colors{i},'FaceAlpha',0.5);
    ylabel('log(count)');
    xlabel('effective dimensionality');
end
set(gca,'XLim',[0 size(d,2)]);
legend(lhs,d_label);
title('DN screen (controls)');


