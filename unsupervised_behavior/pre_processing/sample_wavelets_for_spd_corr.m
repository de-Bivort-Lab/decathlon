
unsupervised_cen = struct('data',[],'strain',[],'IDs',[],'labels',[]);
for i=1:numel(all_cen)
   unsupervised_cen.data =  [unsupervised_cen.data; all_cen{i}.data];
   unsupervised_cen.strain =  [unsupervised_cen.strain; all_cen{i}.strain];
   unsupervised_cen.IDs =  [unsupervised_cen.IDs; all_cen{i}.IDs];
end

%%

unsupervised_cen.labels = cellfun(@(s,i) sprintf('%s_%03.0f',s,i),...
    unsupervised_cen.strain,unsupervised_cen.IDs,'UniformOutput',false);

% 
u_labs = unique(unsupervised_cen.labels);
u_idx = str_list_contains(u_labs,unsupervised_cen.labels);
u_idx = cellfun(@(ui) find(ui,1), num2cell(u_idx,2));

unsupervised_cen.data = unsupervised_cen.data(u_idx);
unsupervised_cen.strain = unsupervised_cen.strain(u_idx);
unsupervised_cen.IDs = unsupervised_cen.IDs(u_idx);
unsupervised_cen.labels = unsupervised_cen.labels(u_idx);

%%

fdir = autoDir;
fpaths = recursiveSearch(fdir,'keyword','projection');
labels = cell(numel(fpaths),1);
parameters = setRunParameters([]);
all_speed = cell(numel(fpaths),1);
all_wavelets = cell(numel(fpaths),1);
npts = 1000;


for i=1:numel(fpaths)
    
    [~,fLabel] = fileparts(fpaths{i});
    j = find(fLabel=='_');
    u_id = sprintf('%s_%03.0f',fLabel(j(end-1)+1:j(end)-1),str2double(fLabel(j(end)+1:end)));
    cen_idx = find(strcmp(unsupervised_cen.labels,u_id));
    if any(cen_idx)
        tmp_cen = unsupervised_cen.data{cen_idx};
        mask = find(mod(1:size(tmp_cen,1),10)==0);
        dec_cen = tmp_cen(mask,:);
        tmp_spd = sqrt(sum([0 0;diff(dec_cen)].^2,2));
        spd = interp1(mask,tmp_spd,1:size(tmp_cen,1))';
        spd(isnan(spd))=0;

        fprintf('loading projections file %i of %i\n',i,numel(fpaths));
        load(fpaths{i},'projections');  

        fprintf('calculating wavelets for file %i of %i\n',i,numel(fpaths));
        [amplitudes,~] = findWavelets(projections,50,parameters);
        wavelet_data = bsxfun(@rdivide,amplitudes,sum(amplitudes,2));

        sample_idx = randperm(size(wavelet_data,1),npts);
        all_speed{i} = spd(sample_idx);
        all_wavelets{i} = wavelet_data(sample_idx,:);

        clear spd tmp_cen tmp_spd wavelet_data projections amplitudes
    end
end