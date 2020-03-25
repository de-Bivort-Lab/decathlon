%% decathlon connected components dimensionality analysis

D = load_decathlon_structs('D13_als');
opts = {'CollapseMode','PCA','CollapseFields','all','Standardize',false};
D = pair_decathlon_structs(D,opts{:});
D = standardize_by_field(D);

d_label = {'D12';'D3'};
figure; hold on;
lhs = gobjects(numel(D),1);
bins = cell(numel(D),1);
colors = {[0 0 1];[1 0.5 0]};
bh = cell(numel(D),1);
for i=1:numel(D)
    d = D(i).data;
    [out,intrinsic_dim,bins{i}]=decathlonConnCompSweep(d,5000);
    cts = histc(intrinsic_dim,1:size(d,2));
    bh{i} = bar(log10(cts),'EdgeColor','none','FaceColor',colors{i},'FaceAlpha',0.5);
end
legend(cat(1,bh{:}),{'D12';'D3'},'Location','NorthWest');
ylabel('log(count)');
xlabel('effective dimensions');



%%
% compute the median bin
peaks = [22,31];
for i=1:numel(D)
    nbins = cellfun(@(b) numel(unique(b)), bins{i});
    peak_bins = find(nbins==peaks(i),1);
    peak_bins = bins{i}{peak_bins};
    f_text = '';
    for j=1:peaks(i)
        f = D(i).fields(peak_bins==j);
        header = sprintf('COMPONENT #%i\n--------------\n',j);
        f = cellfun(@(ff) sprintf('%s\n',ff), f, 'UniformOutput', false);
        f = cat(2,f{:});
        f_text = cat(2,f_text,header,f,sprintf('\n\n'));
    end
    fid = fopen(sprintf('%s/conn_comp_features_%s.txt',pwd,d_label{i}),'W+');
    fwrite(fid,f_text,'char');
    fclose(fid);
end




%%


figure;
for i=1:numel(D)
    subplot(1,numel(D),i);
    dd = D(i).data(randperm(size(D(i).data,1),192),:);
    out=decathlonConnCompDropN(dd);
    imagesc(out);
    axis('equal','tight');
    ylabel('features dropped');
    xlabel('effective dimensionality');
    title(sprintf('%s - ConnComp heatmap',d_label{i}));
    colormap(logjet_cmap);
end

%% unsupervised dimensionality

figure;
d_label = {'D2';'D3'};
for i=1:numel(D_us)
    subplot(1,numel(D_us),i);
    dd = D_us(i).pdfs;
    out=decathlonConnCompDropN(dd);
    imagesc(out);
    axis('equal','tight');
    ylabel('features dropped');
    xlabel('effective dimensionality');
    title(sprintf('unsupervised %s - ConnComp heatmap',d_label{i}));
    colormap(logjet_cmap);
    colorbar;
end

%%

D = D_us;
d_label = {'D12';'D3'};
figure; hold on;
lhs = gobjects(numel(D),1);
bins = cell(numel(D),1);
colors = {[0 0 1];[1 0.5 0]};
bh = cell(numel(D),1);
for i=1:numel(D)
    d = D(i).pdfs;
    [out,intrinsic_dim,bins{i}]=decathlonConnCompSweep(d,1000);
    cts = histc(intrinsic_dim,1:size(d,2));
    bh{i} = bar(log10(cts),'EdgeColor','none','FaceColor',colors{i},'FaceAlpha',0.5);
end
legend(cat(1,bh{:}),{'D12';'D3'},'Location','NorthWest');
ylabel('log(count)');
xlabel('effective dimensions');