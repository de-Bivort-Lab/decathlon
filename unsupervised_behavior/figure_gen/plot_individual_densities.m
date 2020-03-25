function plot_individual_densities(embedding,sigma,numPoints,rangeVals)
% plot density maps for all individual flies separately


% compute densities for each fly
L = numel(embedding.z_data);
densities = zeros(numPoints,numPoints,L);
for i=1:L
    [~,densities(:,:,i)] = findPointDensity(embedding.z_data{i},sigma,numPoints,rangeVals);
end

% plot individual density maps
for i=1:L
    if mod(i-1,5*5)==0
        figure;
    end
    subplot(5,5,mod(i-1,5*5)+1)
    plot_density(densities(:,:,i),ones(size(densities(:,:,i))));
    colorbar('delete');
    title(['Data Set #' num2str(i)],'fontsize',12,'fontweight','bold');
    drawnow;
end