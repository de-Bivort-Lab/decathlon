function options = parse_processing_options(varargin)
% Set default processing options for decathlon data struct functions and append
% to input options if no options are specified. Enforces minimum processing
% options.

% parse inputs
collapse_mode = 'PCA';
collapse_fields = 'none';
do_trim = true;
pcs = 2;
impute_mode = 'als';
for i=1:numel(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'CollapseFields'
                i = i+1;
                collapse_fields = varargin{i};
            case 'CollapseMode'
                i = i+1;
                collapse_mode = varargin{i};
            case 'Trim'
                i = i+1;
                do_trim = varargin{i};
            case 'PCs'
                i = i+1;
                pcs = varargin{i};
            case 'ImputeMode'
                i = i+1;
                impute_mode = varargin{i};
        end
    end
end

options = {'CollapseFields';collapse_fields;'CollapseMode';collapse_mode;...
    'PCs';pcs;'Trim';do_trim;'ImputeMode';impute_mode};