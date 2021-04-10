function cluster_perm = plot_all_pdf_corr(pdfs)

% plot the combined correlation matrix for all flies
fh1=figure;
ah1 = subplot(2,2,1);
trim_feat = any(isnan(pdfs),2);
is_bk = false(size(trim_feat,1),1);
is_bk(1:167) = true;
[~,~,~,zp] = plotCorr(pdfs(~trim_feat,:),...
    'Patch',false,'Options',{},'Parent',ah1);
axis(ah1,'equal','tight');
title(ah1,'Unsupervised PDFs (all flies)');
cluster_perm = zp;

% plot the correlation Bk-iso flies
figure(fh1);
ah2 = subplot(2,2,2);
[~,r1] = plotCorr(pdfs(~trim_feat&is_bk,zp),...
    'Patch',false,'Options',{},'Cluster',false,'Parent',ah2);
title(ah2,'Bk-iso');
axis(ah2,'equal','tight');

% plot the correlation Nex flies
figure(fh1);
ah3 = subplot(2,2,3);
[~,r2] = plotCorr(pdfs(~trim_feat&~is_bk,zp),...
    'Patch',false,'Options',{},'Cluster',false,'Parent',ah3);
title(ah3,'Nex');
axis(ah3,'equal','tight');

% plot the difference between correlation matrices
figure(fh1);
ah4 = subplot(2,2,4);
imagesc(r2-r1,'Parent',ah4);
title(ah4,'\Delta r-values');
colormap(nanticoke);
caxis([-2 2]);
colorbar;
axis('equal','tight');