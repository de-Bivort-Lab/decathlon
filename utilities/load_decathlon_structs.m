function D = load_decathlon_structs(varargin)

% define file path
fpath = ['C:\Users\winsl0w\Documents\decathlon\decathlon_analysis\matrices\'...
    'decathlon_paper\decathlon_final\decathlon_final_data.mat'];
var_name = 'D13_als_med';
if ~isempty(varargin)
    var_name = varargin{1};
end

D = load(fpath,var_name);
fn = fieldnames(D);
D = D.(fn{1});

