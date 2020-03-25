function D = build_DecathlonDataStruct(fDir, varargin)

warning('off');
% parse inputs
data_field_fun = @getDataFields;
culling_permutation = [];
for i=1:numel(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'FieldFunction'
                i=i+1;
                data_field_fun = varargin{i};
            case 'CullingPermutation'
                i=i+1;
                culling_permutation = varargin{i};
        end
    end
end

% search recursively for .mat files
fPaths = recursiveSearch(fDir,'ext','.mat');

% intialize master struct for data
dec = repmat(struct('data',[],'fields',[],'name',[],'day',[],'ID',[],'meta',[]),11,1);
circ = repmat(struct('data',[],'fields',[],'name',[],'day',[],'ID',[],'meta',[]),11,1);

% read in files sequentially
hwb = waitbar(0,'loading files');

for j = 1:length(fPaths)

    hwb = waitbar(j/numel(fPaths),hwb,sprintf('loading file %i of %i',j,numel(fPaths)));
    load(fPaths{j});                    % read in expmt struct
    
    if isprop(expmt,'meta')
        name = expmt.meta.name;                  % query expmt name
        labels_table = expmt.meta.labels_table;
        date = expmt.meta.date;
    else
        name = expmt.Name;
        if strcmpi(name,'Olfaction')
            labels_table = [];
            labels_table.ID = expmt.ID';
            labels_table.Day = expmt.Day;
            if size(labels_table.ID,2) > size(labels_table.ID,1)
                labels_table.ID = labels_table.ID';
            end
        else
            labels_table = expmt.labels_table;
        end
        date = expmt.date;
    end
    
    switch name
        case 'Circadian'
            
            day =  labels_table.Day(1);   % query testing day
            circ(day).ID = [circ(day).ID; labels_table.ID];
            
            % store values in decathlon data struct
            circ(day).name = name;
            circ(day).day = day;

            % extract experiment metrics
            [data,field_names] = feval(data_field_fun,expmt);
            circ(day).fields = field_names;

            % append metrics to values in decathlon data struct
            if isempty(circ(day).data)
                circ(day).data = data;
            else
                fn = fieldnames(data);
                for i = 1:length(fn)
                    circ(day).data.(fn{i}) = [circ(day).data.(fn{i}); data.(fn{i})];
                end
            end
            
            % get meta data for experiment
            meta = assignMetaData(labels_table, name, date);
            if isempty(circ(day).meta)
                circ(day).meta = meta;
            else
                fn = fieldnames(circ(day).meta);
                for i=1:length(fn)
                    circ(day).meta.(fn{i}) = [circ(day).meta.(fn{i}); meta.(fn{i})];
                end
            end

        otherwise
            
            % query testing day
            day =  labels_table.Day(1);   

            if day==1 && strcmpi(name,'Y-Maze')
                name = 'Culling';
            end
            if strcmpi(name,'Basic Tracking') || strcmpi(name,'Arena Circling') 
               name = 'Arena';
                if isprop(expmt,'meta')
                    expmt.meta.name = name;
                else
                    expmt.Name = name;
                end
            end
            
            % store values in decathlon data struct
            dec(day).ID = [dec(day).ID; labels_table.ID];
            dec(day).day = day;
            dec(day).name = name;

            % extract experiment metrics
            [data,field_names] = feval(data_field_fun,expmt);
            dec(day).fields = field_names;

            % append metrics to values in decathlon data struct
            if isempty(dec(day).data)
                dec(day).data = data;
            else
                fn = fieldnames(data);
                for i = 1:length(fn)
                    dec(day).data.(fn{i}) = [dec(day).data.(fn{i}); data.(fn{i})];
                end
            end
            
            % get meta data and append
            meta = assignMetaData(labels_table, name, date);
            if isempty(dec(day).meta)
                dec(day).meta = meta;
            else
                fn = fieldnames(dec(day).meta);
                for i=1:length(fn)
                    dec(day).meta.(fn{i}) = [dec(day).meta.(fn{i}); meta.(fn{i})];
                end
            end
    end
end

delete(hwb);

% sort pre-culling IDs to match post-culling IDs
if ~isempty(culling_permutation)
    f = fieldnames(dec(1).data);
    for i=1:numel(f)
        dec(1).data.(f{i}) = dec(1).data.(f{i})(culling_permutation.culling);
        dec(1).ID = culling_permutation.decathlon;
    end

    f = fieldnames(dec(1).meta);
    for i=1:numel(f)
        dec(1).meta.(f{i}) = dec(1).meta.(f{i})(culling_permutation.culling);
    end
end

% query max number of data points
nFields = 0;
for i = 1:length(dec)
    dec(i).n = length(dec(i).ID);
    nFields = nFields + length(dec(i).fields);
end
cnFields=0;
for i = 1:length(circ)
    circ(i).n = length(circ(i).ID);
    cnFields = cnFields + length(circ(i).fields);
end

% model nuisance variable effects and replaced with residuals if necessary
%[dec,declm] = modelEffects(dec,{'Plate';'Box';'Tray'});
%[circ,circlm] = modelEffects(circ,{'Plate'});

warning('on');

% initialize data mat (observations x variables)
nFields = nFields + cnFields;
n={dec(:).ID};
n=cat(1,n{:});
dMat = NaN(max(n),nFields);
m = cell(size(dMat));
meta = [];
meta.plate = m;
meta.box = m;
meta.tray = m;
meta.day = dMat;
meta.timeofday = dMat;
dFields = cell(nFields,1);
fct = 0;
nDays = max([dec.day circ.day]);

% sort metrics into data matrix
for i=1:nDays
    
    if i<=length(circ) && ~isempty(circ(i).fields)  
        f = circ(i).fields;
        for j=1:length(f)
            fct=fct+1;
            dMat(circ(i).ID,fct) = circ(i).data.(f{j});
            dMat(circ(i).ID(~circ(i).data.filter),fct) = NaN;
            dFields(fct) = {sprintf('%s %s (%i)',circ(i).name,f{j},i)};
            meta = assign_meta_data(meta, circ(i), fct);
        end 
    end
    
    if i<=length(dec) && ~isempty(dec(i).fields)
        f = dec(i).fields;
        for j=1:numel(f)
            fct=fct+1;
            dMat(dec(i).ID,fct) = dec(i).data.(f{j});
            dMat(dec(i).ID(~dec(i).data.filter),fct) = NaN;
            dFields(fct) = {sprintf('%s %s (%i)',dec(i).name,f{j},i)};
            meta = assign_meta_data(meta, dec(i), fct);
        end 
    end
end

% remove completely empty observations
empty_rows = all(isnan(dMat),2);
dMat(empty_rows,:) = [];
f = fieldnames(meta);
for i=1:numel(f)
    meta.(f{i})(empty_rows,:) = [];
end
n = get_pairwise_sampling(dMat);

% initialize Decathlon data struct
D.data = dMat;
D.fields = dFields;
D.meta = meta;
D.n = n;

% assign meta data from dec struct to meta data matrix
function meta = assign_meta_data(meta, dec, fct)

requested_fields = sort(fieldnames(meta));
existing_fields = sort(fieldnames(dec.meta));
fq = existing_fields(any(str_list_contains(lower(existing_fields), requested_fields),2));
fa = requested_fields(cellfun(@(rf) any(strcmpi(fq,rf)), requested_fields));

for i=1:numel(fa)
    meta.(fa{i})(dec.ID, fct) = dec.meta.(fq{i});
end