
cluster_ranges = [{1268:1341};{1342:1989};{2082:2121};{2286:2532}];
cts = NaN(numel(cluster_ranges,4));
all_modes_present = cell(numel(cluster_ranges),1);
figure;
for i=1:numel(cluster_ranges)
    
    cluster_labels = all_labels(cluster_ranges{i});
    a=regexp(cluster_labels,'P\(mode_{[0-9]*})','match');
    a = cat(1,a{:});
    modes_present = regexp(a,'[0-9]*','match');
    modes_present = cellfun(@str2double,modes_present);
    all_modes_present{i} = modes_present;
    
    subplot(2,2,i);
    imagesc(pdfs(:,modes_present)');
    molaspass=interp1([1 51 102 153 204 256],[0 0 0; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 1 1 1],1:256);
    colormap(molaspass);
    caxis([0 1]);
    colorbar;

    a=regexp(cluster_labels,'P\([0-9]*,[0-9]*','match');
    a = cat(1,a{:});
    from = regexp(a,'[0-9]*(?=,)','match');
    from = cellfun(@str2double,cat(1,from{:}));
    to = regexp(a,'(?<=,)[0-9]*','match');
    to = cellfun(@str2double,cat(1,to{:}));
    from_cts = histc(from,modeIDs);
    from_cts = from_cts./sum(from_cts);
    to_cts = histc(to,modeIDs);
    to_cts = to_cts./sum(to_cts);
%     subplot(2,2,i);
%     bar([from_cts to_cts],'stacked');
    from_ct = ismember(from,modes_present);
    to_ct = ismember(to,modes_present);
    p = numel(from_ct)*((1/numel(modeIDs))*numel(modes_present));
    q = numel(from_ct)*((1/numel(modeIDs))*(numel(modeIDs)-numel(modes_present)));
    p2p = p*(1/numel(modeIDs))*(numel(modes_present)-1);
    q2q = q*((1/numel(modeIDs))*(numel(modeIDs)-numel(modes_present)-1));
    
    
    cts(i,1) = sum(from_ct)/p;
    cts(i,2) = sum(~from_ct)/q;
    cts(i,3) = sum(to_ct)/p;
    cts(i,4) = sum(~to_ct)/q;
    cts(i,5) = sum(from_ct & to_ct)/p2p;
    cts(i,6) = sum(~from_ct & ~to_ct)/q2q;
end


%
figure;
bh=bar(log(cts));
leg_labels = {'transitions A-x';'transitions B-x';'transitions x-A';...
    'transitions x-B';'transitions A-A';'transitions B-B'};
legend(bh,leg_labels);
ylabel('log(obs/exp)');
xlabel('cluster #');

%
modes_present = cat(1,all_modes_present{:});
modes_not_present = modeIDs(~ismember(modeIDs,modes_present));
figure;
imagesc(pdfs(:,modes_not_present)');
molaspass=interp1([1 51 102 153 204 256],[0 0 0; 0 0 .75; .5 0 .8; 1 .1 0; 1 .9 0; 1 1 1],1:256);
colormap(molaspass);
colorbar;