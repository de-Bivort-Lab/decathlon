function fig_handle = plotCorr_byDate(dMat,dFields,varargin)
% plot decathlon correlation matrix sorted by day of testing

% get unique day numbers
[~,~,days] = parse_fieldnames(dFields);
udays = unique(days);

% sort by data and field labels by day of testing
sorted_dMat = [];
sorted_dFields = [];

for i=1:numel(udays) 
    sorted_dMat = [sorted_dMat dMat(:,days==udays(i))];
    sorted_dFields = [sorted_dFields; dFields(days==udays(i))];
end

[fig_handle,~,~] = plotCorr(sorted_dMat,...
    'Labels',sorted_dFields,'Cluster',false,'Signed',false,varargin{:});

