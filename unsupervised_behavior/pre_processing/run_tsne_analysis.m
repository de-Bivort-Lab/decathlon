
% select directory containing .mat files containing embedding data
[fname,fpath] = uigetfile('Select decathlon_unsupervised_embedding_timeseries.mat');
embedding = load(cat(2,fpath,fname));
embedding = embedding.embeddings;

% find t-SNE speed threshold for parsing t-SNE trajectories into pauses at modes
make_plot = true;
[sigma,embedding.z_thresh] = fit_tsne_z_logspeed_gmm(embedding.z_speed,make_plot);

% set kernel parameters for density estimation
L = numel(embedding.z_data);
maxVal = max(max(abs(cat(1,embedding.z_data{:}))));
maxVal = round(maxVal * 1.1);
sig = sigma(1);
numPoints = 501;
rangeVals = [-maxVal maxVal];

% compute raw density without any filtering
[~,density] = findPointDensity(cat(1,embedding.z_data{:}),sig,numPoints,rangeVals);
figure('Name','Combined density maps');
subplot(1,2,1);
plot_density(density,ones(size(density)));
title('Unfiltered density map');

% filter data for pts pausing below the speed threshold for a min duration
frame_rate = 100;
min_dur = 0.05;
slow_mode = cellfun(@(zz) zz<embedding.z_thresh, embedding.z_speed, 'UniformOutput', false);
filt_slow_mode = cellfun(@(sm) filter_speed_mask(sm,frame_rate*min_dur),slow_mode, 'UniformOutput', false);
filt_embedded_pts = cellfun(@(emv,fsm) emv(fsm,:), embedding.z_data, filt_slow_mode,'UniformOutput',false);

% generate density from filtered data
[xx,filt_density] = findPointDensity(cat(1,filt_embedded_pts{:}),sig,numPoints,rangeVals);
subplot(1,2,2);
plot_density(filt_density,ones(size(density)));
title('Slow speed filtered density map');

% uncomment to plot individual densities
% plot_individual_densities(embedding,sig,numPoints,rangeVals);

% uncomment to plot densities by genotype and t-SNE speed mode
% plot_genotype_densities(embedding,sig,numPoints,rangeVals,[0 max(density(:))*1.02]);

% uncomment to plot individual z-trace position samples
% plot_tsne_position_samples(embedding);


%% Assign behavioral mode classifications

% watershed density map and classify each frame as one of the resulting modes
[frameIDs,modeIDs,pdfs,idxMap] = mode_from_embeddingValues(filt_density,xx,embedding.z_data);

% plot density for BK-iso
is_bk = false(size(frameIDs,1),1);
is_bk(1:166) = true;
slow_embedded_pts = cellfun(@(emv,fsm) emv(fsm,:), embedding.z_data, slow_mode,'UniformOutput',false);
bk_slow_pts = cat(1,slow_embedded_pts{is_bk});
[xx,bk_density] = findPointDensity(bk_slow_pts,sig,numPoints,rangeVals);
figure;
plot_density(bk_density,idxMap,'Numbered',0,'OutlineDensity',filt_density);

% plot individual pdfs and mode map pdf
plot_mode_pdfs(idxMap,pdfs);
figure;
plot_density(filt_density,idxMap,'Numbered',1);

plot_all_pdf_corr(pdfs);

% read in labels for each mode from hand annotated file
mode_label_path = 'D:\D2D3_combined_unsupervised\master data\mode_descriptions.txt';
mode_labels = load_manual_mode_labels(mode_label_path);

% plot PDF correlation matrices for Nex and Bk-iso separately
[fh,~,~,zp]=plotCorr(pdfs(is_bk,:),'Labels',mode_labels(:),'Patch',false,'FontSize',8);
close(fh(2));
title('Bk-iso PDF corrmat');
axis('equal','tight');
fh = plotCorr(pdfs(~is_bk,zp),'Labels',mode_labels(zp),'Patch',false,'Cluster',false,'FontSize',8);
close(fh(2));
title('Nex PDF corrmat');
axis('equal','tight');

%% use classifications to assign pause bouts in each mode (over some min duration)

% filter frame labels by slow bout mask
allModes = unique(idxMap);
targetDuration = 0.05*100;
filt_frameIDs = frameIDs;
for i=1:numel(filt_slow_mode)
   filt_frameIDs{i}(~filt_slow_mode{i}) = 0; 
end

% parse frame classifications into bouts of pauses at mode peaks
nflies = numel(filt_frameIDs);
[starts,stops,durations,sampleFrames,sampleDurations] = cellfun(@(x) ...
    modeBouts(x,modeIDs,frame_rate*min_dur),filt_frameIDs,'UniformOutput',false);

%% use mode bout assignments to calculate transition probabilities between modes

% calculate the transition probabilities between modes
mode_trans_probs = cellfun(@(si,i) find_mode_transitions(si,modeIDs,i,numel(starts)),...
    starts,num2cell(1:numel(starts))','UniformOutput',false);

% compute average transition prob matrix from all flies
all_trans = cat(3,mode_trans_probs{:});
avg_trans = nanmean(all_trans(:,:,is_bk),3);
figure; imagesc(avg_trans(zp,zp));
ylabel('initial state');
xlabel('final state');
title('Unsupervised transition probabilities');
axis('equal','tight')
caxis([0 .15]);
colormap(density_cmap);
colorbar;

flat_trans_probs = cellfun(@(mtp) mtp(:), mode_trans_probs,'UniformOutput',false);

all_features = [pdfs cat(2,flat_trans_probs{:})'];
%mode_labels = arrayfun(@(i) sprintf('P(mode_{%i})',i), modeIDs, 'UniformOutput', false)';
ii = repmat(modeIDs,numel(modeIDs),1);
jj = repmat(modeIDs',1,numel(modeIDs));
trans_labels = arrayfun(@(i,j) sprintf('P(%i,%i)',i,j), ii(:), jj(:), ...
    'UniformOutput',false);

trim_feat = any(isnan(all_features),2);
[~,~,~,zp] = plotCorr(all_features(~trim_feat,:),'Patch',false,'Options',{});
title('Unsupervised features correlation matrix (all flies)');
all_labels = [mode_labels; trans_labels];
all_labels = all_labels(zp);

fh1=figure;
ah1 = subplot(1,2,1);
plotCorr(all_features(~trim_feat&is_bk,zp),...
    'Patch',false,'Options',{},'Cluster',false,'Parent',ah1);
title(ah1,'Unsupervised features correlation matrix (Bk-iso-1)');
axis(ah1,'equal','tight');
figure(fh1);
ah2 = subplot(1,2,2);
plotCorr(all_features(~trim_feat&~is_bk,zp),...
    'Patch',false,'Options',{},'Cluster',false,'Parent',ah2);
title(ah2,'Unsupervised features correlation matrix (Nex)');
axis(ah2,'equal','tight');

%% format transition probabilities and PDFs into decathlon structs

labels = regexp(embedding.label,'_','split');
strains = cellfun(@(s) s{1}, labels, 'UniformOutput', false);
ids = cellfun(@(s) str2double(s{2}), labels);
D_us = repmat(struct('strain',[],'ID',[],'pdfs',[],'pdf_labels',[]...
    ,'trans',[],'trans_labels',[]),2,1);
for i=1:2
   if i<2
       mask = is_bk;
   else
       mask = ~is_bk;
   end
   D_us(i).strain = strains(mask); 
   D_us(i).ID = ids(mask); 
   D_us(i).pdfs = pdfs(mask,:);
   D_us(i).pdf_labels = mode_labels; 
   D_us(i).trans = cat(2,flat_trans_probs{mask})';
   D_us(i).trans_labels = trans_labels; 
end

% store a condensed version of the matrix without transition probs
D_u = repmat(struct('data',[],'fields',[],'ID',[],'strain',[]),2,1);
for i=1:numel(D_us)
    D_u(i).data = D_us(i).pdfs;
    D_u(i).fields = mode_labels; 
    D_u(i).ID = D_us(i).ID;
    D_u(i).strain = D_us(i).strain;
end


%%

% plot the pdf of occupancy durations
durations = cat(2,durations{:});
durations = num2cell(durations,2);
plot_mode_occupancy_pdf(durations);

% create sample movies of each mode
sampleFrames =cat(2,sampleFrames{:});
sampleFrames = num2cell(sampleFrames,2);
sampleDurations = cat(2,sampleDurations{:});
sampleDurations = num2cell(sampleDurations,2);

% generate sampling vector to generate tiled movie for each mode
movieVec = cellfun(@getModeMovieVector,sampleFrames,sampleDurations,'UniformOutput',false);

% make tiled movies
aligned_movie_dir = {'E:\decathlon unsupervised videos\decathlon 2-2018 aligned video';...
    'E:\decathlon unsupervised videos\decathlon 1-2019 aligned';'I:\decathlon 1-2019 aligned video'};
make_all_tiled_mode_movies(movieVec,aligned_movie_dir,numel(frameIDs));
