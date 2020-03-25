function plot_density_by_genotype(embedding,sig,numPoints,rangeVals,clim)

% define unique strains
strains = {'Bk-iso-1';'Nex'};
figure('Name','Genetic Background Density Maps');
split_densities = cell(3,numel(unique(embedding.strain)));
for i=1:numel(strains)
    ah = subplot(3,numel(strains),i);
    
    % get density estimate from embedding data
    f = str_list_contains(embedding.label,strains{i});
    embedded_pts = combineCells(embedding.z_data(f));
    embedded_z = cat(1,embedding.z_speed{f});
    [xx,split_densities{1,i}] = findPointDensity(embedded_pts,sig,numPoints,rangeVals);

    % plot total density map
    plot_density(split_densities{1,i},ones(numPoints),'CLim',clim);
    title(sprintf('%s (unfiltered)',strains{i}));
    
    % get density estimate from embedding data
    ah = subplot(3,numel(strains),i+2);
    ff = embedded_z < embedding.z_thresh;
    [xx,split_densities{2,i}] = findPointDensity(embedded_pts(ff,:),sig,numPoints,rangeVals);

    % plot total density map
    plot_density(split_densities{2,i},ones(numPoints),'CLim',clim);
    title(sprintf('%s (slow mode)',strains{i}));

    
    % get density estimate from embedding data
    ah = subplot(3,numel(strains),i+4);
    ff = embedded_z > embedding.z_thresh;
    [xx,split_densities{3,i}] = findPointDensity(embedded_pts(ff,:),sig,numPoints,rangeVals);

    % plot total density map
    plot_density(split_densities{3,i},ones(numPoints),'CLim',clim);
    title(sprintf('%s (fast mode)',strains{i}));
end