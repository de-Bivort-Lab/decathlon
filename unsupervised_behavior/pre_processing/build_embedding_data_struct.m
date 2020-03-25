function embedding = build_embedding_data_struct(fDir)
% load t-SNE embedding data from .mat files and use file names to sort
% data by genotype and fly ID number
    
% search for and load all embedding files
fPaths = recursiveSearch(fDir,'keyword','embedding','ext','.mat');
embeddingValues = cell(numel(fPaths),1);
for i =1:length(fPaths)
    fprintf('loading embeddings file %i of %i\n',i,numel(fPaths));
    load(fPaths{i},'emVal');
    embeddingValues(i) = {emVal};
end

% parse genotype and fly ID from filenames
[~,fLabels,~] = cellfun(@fileparts,fPaths,'UniformOutput',false);
strain_idx = cellfun(@(s) find(s=='_'),fLabels,'UniformOutput',false);
strain_idx = cat(1,strain_idx{:});
id_idx = strain_idx(:,end);
strain_idx = strain_idx(:,end-1);
unique_strain = cellfun(@(s,ii,jj) s(ii+1:jj-1), ...
    fLabels, num2cell(strain_idx), num2cell(id_idx), 'UniformOutput', false);
unique_IDs = cellfun(@(s,ii) str2double(s(ii+1:end)), ...
    fLabels, num2cell(id_idx));

% generate a unique identying label for each file and find sorting permutation
unique_label = cellfun(@(s,ii) sprintf('%s_%03.0f',s,ii),...
    unique_strain,num2cell(unique_IDs),'UniformOutput',false);
[~,p] = sort(unique_label);

% initialze embedding data struct
embedding = struct('z_data',[],'z_speed',[],'mode_data',[],'mode_labels',[],...
    'mode_pdfs',[],'mode_bouts',[],'label',[],'strain',[],'ID',[]);

% permute labels and embedding data
embedding.z_data = embeddingValues(p);
embedding.label = unique_label(p);
embedding.strain = unique_strain(p);
embedding.ID = unique_IDs(p);

% compute embedding trajectory speed
embedding.z_speed = cellfun(@(x) sqrt([0;diff(x(:,1))].^2 + [0;diff(x(:,2))].^2).*100,...
    embedding.z_data,'UniformOutput',false);