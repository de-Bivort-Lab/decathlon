%% Get parent directory of all decathlon files

[fDir] = uigetdir('C:\Users\debivort\Documents\MATLAB\Decathlon Raw Data',...
    'Select directory containing expmt structs to be analyzed');

% intialize master struct for data
dec = repmat(struct('data',[],'fields',[],'name',[],'day',[],'ID',[]),11,1);
circ = repmat(struct('data',[],'fields',[],'name',[],'day',[],'ID',[]),11,1);

fPaths = recursiveSearch(fDir);

%% read in files sequentially

hwb = waitbar(0,'loading files');

for j = 1:length(fPaths)
    
    hwb = waitbar(j/length(fPaths),hwb,['loading file ' num2str(j) ' of ' num2str(length(fPaths))]);
    load(fPaths{j});                    % read in expmt struct
    name = expmt.Name;                  % query expmt name
    
    switch name
        case 'Circadian'
            
            day =  expmt.labels_table.Day(1);   % query testing day
            circ(day).ID = [circ(day).ID; expmt.labels_table.ID];
            expmt.nTracks = length(expmt.labels_table.ID);
            
            % store values in decathlon data struct
            circ(day).name = name;
            circ(day).day = day;

            % extract experiment metrics
            [data,field_names] = getDataFields_legacy(expmt);
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
            
        otherwise
            
            switch name
                case 'Olfaction'
                    day = expmt.Day;
                    dec(day).ID = [dec(day).ID; expmt.ID];
                otherwise
                    day =  expmt.labels_table.Day(1);   % query testing day
                    dec(day).ID = [dec(day).ID; expmt.labels_table.ID];
                    expmt.nTracks = length(expmt.labels_table.ID);
            end

            % store values in decathlon data struct
            dec(day).name = name;
            dec(day).day = day;

            % extract experiment metrics
            [data,field_names] = getDataFields_legacy(expmt);
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
            
    end
    
    
end
    
delete(hwb);

%% create data matrix and create labels array

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
%[dec,declm] = modelEffects(dec,nFields,'TimeofDay',false);
%[circ,circlm] = modelEffects(circ,cnFields,'TimeofDay',false);

% initialize data mat (observations x variables)
nFields = nFields + cnFields;
n={dec(:).ID};
n=cat(1,n{:});
dMat = NaN(max(n),nFields);
dFields = cell(nFields,1);
fct = 0;
nDays = max([dec.day circ.day]);

for i=1:nDays
    
    if i<=length(circ) && ~isempty(circ(i).fields)
        
        f = circ(i).fields;

        for j=1:length(f)

            fct=fct+1;
            dMat(circ(i).ID,fct) = circ(i).data.(f{j});
            dMat(circ(i).ID(~circ(i).data.filter),fct) = NaN;
            dFields(fct) = {[circ(i).name ' ' f{j} ' (' num2str(i) ')']};
            
            if any(dMat(:,end)==0)
                disp(i);
            end
        end 
    end
    
    if i<=length(dec) && ~isempty(dec(i).fields)
        
        f = dec(i).fields;

        for j=1:length(f)

            fct=fct+1;
            dMat(dec(i).ID,fct) = dec(i).data.(f{j});
            dMat(dec(i).ID(~dec(i).data.filter),fct) = NaN;
            dFields(fct) = {[dec(i).name ' ' f{j} ' (' num2str(i) ')']};
        end 
    end
    
    
end

% calculate number of samples for each pairwise comparison
n = NaN(nFields);
for i=1:nFields
    for j=1:nFields       
        if i~=j
            n(i,j) = sum(~isnan(dMat(:,i)) & ~isnan(dMat(:,j)));
        end
    end
end

% delete empty struct entries
dec(arrayfun(@(x) isempty(dec(x).data),1:length(dec)))=[];
circ(arrayfun(@(x) isempty(circ(x).data),1:length(circ)))=[];


dMat(~any(~isnan(dMat')),:)=[];
clearvars -except n dec circ dMat dFields nFields declm circlm
