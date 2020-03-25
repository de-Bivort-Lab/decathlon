

slow_embedded_pts = cellfun(@(emv,fsm) emv(fsm,:), embedding.z_data, slow_mode,'UniformOutput',false);
bk_slow_pts = cat(1,slow_embedded_pts{is_bk});
[xx,bk_density] = findPointDensity(bk_slow_pts,sig,numPoints,rangeVals);

figure;
plot_density(bk_density,ones(size(bk_density)),'Numbered',0,'OutlineDensity',ones(size(bk_density)));
figure;
plot_density(bk_density,idxMap,'Numbered',0,'OutlineDensity',filt_density);


% plot shared pdf
pdf_bk = pdfs(is_bk,:);
plot_mode_pdfs(idxMap,pdf_bk);

% plot some sample modes
idx = [42 84 126];
for i=1:numel(idx)
    figure;
    plot_mode_pdfs(idxMap,pdf_bk(idx(i),:));
end



% plot matrices separately
[fh,~,~,~]=plotCorr(pdfs(is_bk,zp(pp)),'Labels',mode_labels(zp(pp)),'Patch',false,'FontSize',12,'Cluster',false);
close(fh(2));
title('Bk-iso PDF corrmat');
axis('equal','tight');
fh = plotCorr(pdfs(~is_bk,zp(pp)),'Labels',mode_labels(zp(pp)),'Patch',false,'Cluster',false,'FontSize',12);
close(fh(2));
title('Nex PDF corrmat');
axis('equal','tight');