% Compare variances of behavioral and gene expression data sets in D1 and D2
fdir = ['D:\decathlon_preprint_code_data_figures\decathlon_analysis\'...
    'matrices\decathlon_paper\decathlon_final\'];
D = load_decathlon_structs(fdir,'D123_als');
D_p = pair_decathlon_structs(D,'CollapseMode','PCA','CollapseFields','none');
[coef,~,~,~,v1,mu] = pca(D_p(1).data);
x = (D_p(2).data - mu)*coef;
x = x(~all(isnan(x),2),:);
v2 = var(x)';
v2 = (v2./sum(v2)).*100;

% plot variances
figure;
subplot(1,3,1); hold on;
lh1 = plot(v1,'k--','LineWidth',1.5);
lh2 = plot(v2,'r-','LineWidth',1.5);
xlabel('PC no.');
ylabel('% variance explained');
set(gca,'YScale','log');
legend([lh1,lh2],{'D1';'D2'});
title('Full Matrices - behavior');

% generate plot for distilled matrix
D_p = pair_decathlon_structs(D,'CollapseMode','PCA','CollapseFields','all');
[coef,~,~,~,v1,mu] = pca(D_p(1).data);
x = (D_p(2).data - mu)*coef;
x = x(~all(isnan(x),2),:);
v2 = var(x)';
v2 = (v2./sum(v2)).*100;

% plot variances
subplot(1,3,2); hold on;
lh1 = plot(v1,'k--','LineWidth',1.5);
lh2 = plot(v2,'r-','LineWidth',1.5);
xlabel('PC no.');
ylabel('% variance explained');
set(gca,'YScale','log');
legend([lh1,lh2],{'D1';'D2'});
title('Distilled Matrices - behavior');

%% Sequencing data

N = 20;
D_p = pair_rnaseq(D_seq(1:2));
for i=1:numel(D_p)
   D_p(i).data = D_p(i).data(all(~isnan(D_p(i).data),2),:);
   D_p(i).data = D_p(i).data(sum(D_p(i).data,2)>.5E6,:);
   qn = quantile_normalize(D_p(i).data');
   D_p(i).data = qn';
end
[coef,~,~,~,v1,mu] = pca(D_p(1).data,'NumComponents',N);
x = (D_p(2).data - mu)*coef;
x = x(~all(isnan(x),2),:);
v2 = var(x)';
v2 = (v2./sum(v2)).*100;

% plot variances
subplot(1,3,3); hold on;
lh1 = plot(v1(1:N),'k--','LineWidth',1.5);
lh2 = plot(v2(1:N),'r-','LineWidth',1.5);
xlabel('PC no.');
ylabel('% variance explained');
set(gca,'YScale','log');
legend([lh1,lh2],{'D1';'D2'});
title('Gene Expression');