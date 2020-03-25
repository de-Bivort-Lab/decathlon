
% define idle mode numbers
idle_modes = uint8([32 39 43 67 68 69]);
idle_mode_frames = cellfun(@(ffid) find(ismember(ffid,idle_modes)), ...
    filt_frameIDs, 'UniformOutput', false);
idle_mode_frames = cellfun(@(imf) imf(randperm(numel(imf),100)), ...
    idle_mode_frames, 'UniformOutput', false);

%%

% initialize projection file paths
fdir = autoDir;
fpaths = recursiveSearch(fdir,'keyword','projection');
labels = embedding.label;
fLabels = regexp(fpaths,'[_|.]','split');
fLabels = cellfun(@(fl) sprintf('%s_%03.0f',fl{8},str2double(fl{9})), ...
    fLabels, 'UniformOutput', false);
[~,p] = sort(fLabels);
fpaths = fpaths(p);
parameters = setRunParameters([]);
idle_wavelets = cell(numel(fpaths),1);

%%

for i=1:numel(fpaths)
    
    fprintf('loading projections file %i of %i\n',i,numel(fpaths));
    load(fpaths{i},'projections');

    fprintf('calculating wavelets for file %i of %i\n',i,numel(fpaths));
    [amplitudes,~] = findWavelets(projections,50,parameters);
    wavelet_data = bsxfun(@rdivide,amplitudes,sum(amplitudes,2));
    idle_wavelets{i} = wavelet_data(idle_mode_frames{i},:);

    clear wavelet_data projections amplitudes
end

%%

fprintf(1,'Finding t-SNE Embedding for the Training Set\n');
trainingSetData = cat(1,idle_wavelets{:});
[trainingEmbedding,betas,P,errors] = run_tSne(trainingSetData,parameters);

%%

% define idle mode numbers
idle_modes = uint8([32 39 43 67 68 69]);
idle_mode_frames = cellfun(@(ffid) find(ismember(ffid,idle_modes)), ...
    filt_frameIDs, 'UniformOutput', false);
out_path = 'D:\D2D3_combined_unsupervised\D2D3_combined_embeddings\idle_embeddings\';
L = numel(fpaths);

fprintf(1,'Finding t-SNE Embedding for each file\n');
num_points = 15000;
parameters.embedding_batchSize = num_points;
for i=116:numel(fpaths)
    
    fprintf(1,'\t Finding Embbeddings for File #%4i out of %4i\n',i,L);
    load(fpaths{i},'projections');
    [~,fname]=fileparts(fpaths{i});

    fprintf('calculating wavelets for file %i of %i\n',i,numel(fpaths));
    [amplitudes,~] = findWavelets(projections,50,parameters);
    clear projections
    
    wavelet_data = bsxfun(@rdivide,amplitudes,sum(amplitudes,2));
    n = num_points;
    n(n>numel(idle_mode_frames{i})) = numel(idle_mode_frames{i});
    idle_idx = randperm(numel(idle_mode_frames{i}),n);
    idle_wavelets = wavelet_data(idle_mode_frames{i}(idle_idx),:);
    em = findTDistributedProjections_fmin(...
        idle_wavelets,trainingSetData,trainingEmbedding,parameters);
    
    save(sprintf('%s%s.mat',out_path,fname),'em');

    clear wavelet_data amplitudes em idle_wavelets
end

