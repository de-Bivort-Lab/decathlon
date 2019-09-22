% pre-process decathlon persistence data consistent with decathlon
% pre-processing

fDir = autoDir;
if ~iscell(fDir)
   fDir = {fDir}; 
end
fPaths = cellfun(@(d) getHiddenMatDir(d), fDir, 'UniformOutput',false);
fPaths = cat(2,fPaths{:});
fPaths = fPaths';
options = {'Raw',true,'Regress',false,...
    'Slide',false,'Bootstrap',false,'Save',false,'Zip',false};


%%

assays = {'Y-Maze';'LED Y-Maze';'Circadian';'Slow Phototaxis';'Optomotor'};
data = cell(10000,1);
d_fields = data;
data(:) = {NaN(192,1)};
d_fields(:) = {''};

for i=1:numel(fPaths)
    fprintf('processing file %i of %i\n',i,numel(fPaths));
    [tmp_dir,~,~] = fileparts(fPaths{i});
    options = {'Dir',[tmp_dir '\'],'Raw',true,'Regress',false,...
        'Slide',true,'Bootstrap',false,'Save',true,'Zip',false};
    if any(str_list_contains(fPaths(i),{'Y-maze';'LED'}))
        load(fPaths{i});
        tmp_data = parse_v2_ymaze_data(flyTracks);
        tmp_data.nTrials = tmp_data.n;
        tmp_data.hand_clumpiness = tmp_data.clumpiness;
        tmp_data.hand_switchiness = tmp_data.switchiness;
        tmp_data = rmfield(tmp_data,{'n';'clumpiness';'switchiness';'t';'idx'});
        fields = fieldnames(tmp_data);
        
        % get meta data
        labels = flyTracks.labels;
        assay = flyTracks.exp;
    else
        expmt = analyze_multiFile(options{:});
        [tmp_data,fields] = getDataFields_legacy(expmt);
        
        % get meta data
        labels = expmt.labels_table;
        assay = expmt.Name;
    end
    var_names = labels.Properties.VariableNames;
    if any(strcmp(var_names,'Day'))
        day = labels.Day(1);
        ids = labels.ID(~isnan(labels.ID));
    else
        day = labels.day(1);
        ids = labels.ID(~isnan(labels.ID));
    end
        
    for j=1:numel(fields)
       tmp_fn = sprintf('%s %s (%i)',assay,fields{j},day);
       fidx = find(strcmpi(d_fields,tmp_fn),1);
       if isempty(fidx)
           fidx = find(cellfun(@isempty,d_fields),1);
       end
       data{fidx}(ids) = tmp_data.(fields{j})(1:numel(ids));
       d_fields{fidx} = tmp_fn;
    end
    
    close all
end

%%

empty_cells = cellfun(@isempty,d_fields);
d_fields(empty_cells) = [];
data(empty_cells) = [];

%%

