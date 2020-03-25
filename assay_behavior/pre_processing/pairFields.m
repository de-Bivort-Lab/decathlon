function [A_sorted,B_sorted] = pairFields(A,B,varargin)
% Uses the field names in decathlon data structs A and B to remove any fields 
% not common to both structures and sorts the remaining fields to match. 
% Name-Value pair Trim specifies whether to trim Day info in the field
% names prior to sorting (with the execption of Circadian measures which are
% measured on every day of testing. Specifying Trim = true will sort 
% non-circadian measures regardless of day of testing and will sort circadian 
% measures by day of testing.

trim=false;
for i=1:length(varargin)
    
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'Trim'
                i=i+1;
                trim = varargin{i};
        end
    end
end

% get the field names of each dataset
f1 = standardize_fieldnames(A.fields);
f2 = standardize_fieldnames(B.fields);

if trim
    % remove day of testing information from non-circadian measures
    af1 = find(cellfun(@isempty,strfind(f1,'Circadian')));
    trimmed = cellfun(@(x,y) x(1:y-2),f1(af1),strfind(f1(af1),'('),'UniformOutput',false);
    f1(af1(~cellfun(@isempty,trimmed))) = trimmed(~cellfun(@isempty,trimmed));
    af2 = find(cellfun(@isempty,strfind(f2,'Circadian')));
    trimmed = cellfun(@(x,y) x(1:y-2),f2(af2),strfind(f2(af2),'('),'UniformOutput',false);
    f2(af2(~cellfun(@isempty,trimmed))) = trimmed(~cellfun(@isempty,trimmed));
end

% get mapping between datasets A and B
permutation = cellfun(@(x) find(strcmpi(x,f2(:)),1,'Last'),f1,'UniformOutput',false);
idxA = 1:length(f1);
idxA(cellfun(@isempty,permutation))=[];
permutation(cellfun(@isempty,permutation))=[];
permutation = cat(1,permutation{:});

% permute data and fields
A_sorted.data = A.data(:,idxA);
A_sorted.fields = f1(idxA);
A_sorted.n = get_pairwise_sampling(A_sorted.data);
meta_f  = fieldnames(A.meta);
for i=1:numel(meta_f)
    A_sorted.meta.(meta_f{i}) = A.meta.(meta_f{i})(:,idxA);
end
if isfield(A,'imputed')
   A_sorted.imputed = A.imputed(:,idxA); 
end

% sort any remaining fields
fn = fieldnames(A);
fn = fn(~ismember(fn,{'fields';'data';'meta';'imputed'}));
nf_a = numel(A.fields);
for i=1:numel(fn)
   tmp = A.(fn{i});
   f_dim = find(size(tmp)==nf_a);
   if numel(f_dim)==1
       switch f_dim
           case 1
               A_sorted.(fn{i}) = tmp(idxA,:);
           case 2
               A_sorted.(fn{i}) = tmp(:,idxA);
       end
   elseif numel(f_dim)==2
      A_sorted.(fn{i}) = tmp(idxA,idxA);
   end
end

B_sorted.data = B.data(:,permutation);
B_sorted.fields = f2(permutation);
B_sorted.n = get_pairwise_sampling(B_sorted.data);
meta_f  = fieldnames(B.meta);
for i=1:numel(meta_f)
    B_sorted.meta.(meta_f{i}) = B.meta.(meta_f{i})(:,permutation);
end
if isfield(B,'imputed')
   B_sorted.imputed = B.imputed(:,permutation); 
end

% sort any remaining fields
fn = fieldnames(B);
fn = fn(~ismember(fn,{'fields';'data';'meta';'imputed'}));
nf_b = numel(B.fields);
for i=1:numel(fn)
   tmp = B.(fn{i});
   f_dim = find(size(tmp)==nf_b);
   if numel(f_dim)==1
       switch f_dim
           case 1
               B_sorted.(fn{i}) = tmp(permutation,:);
           case 2
               B_sorted.(fn{i}) = tmp(:,permutation);
       end
   elseif numel(f_dim)==2
      B_sorted.(fn{i}) = tmp(permutation,permutation);
   end
   
end

