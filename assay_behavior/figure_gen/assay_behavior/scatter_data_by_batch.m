function scatter_data_by_batch(D)
pd = makedist('Normal','mu',0,'sigma',0.08);
opts = {'Marker'; 'o'; 'LineStyle'; 'none';...
    'MarkerFaceColor'; 'k'; 'MarkerEdgeColor'; 'none';...
    'MarkerSize'; 1.5; 'LineWidth'; 1};


for i=1:numel(D)
    figure('Name',sprintf('Struct %i',i));
    data = D(i).data;
    assays = parse_fieldnames(D(i).fields);
    labels = pretty_labels(D(i).fields);
    
    for j = 1:size(data,2)
       timeofday = D(i).meta.timeofday(:,j);
       metric = data(:,j);
       ax = subplot(9,9,j);
       hold(ax,'on');
       if strcmpi(assays{j},'Olfaction')
           circ_idx = find(strcmpi(assays(1:j),'Circadian'),1,'Last');
           plate = D(i).meta.plate(:,circ_idx);
           plate(cellfun(@isempty,plate)) = {''};
           grps = unique(plate);
           [~,timeofday] = find(str_list_contains(plate,grps));
       end
       grps = unique(timeofday);
       for k=1:numel(grps)
          mask = timeofday == grps(k);
          grp_data = metric(mask);
          xx = repmat(random(pd,1,sum(mask)), 1, 1);
          plot(k+xx',metric(mask),opts{:});
       end
       set(ax,'XTick',[],'YTick',[]);
       title(labels(j),'FontSize',6,'FontWeight','normal');
       hold(ax,'off');
       data(:,j) = metric;
    end
    D(i).data = data;
end