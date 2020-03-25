function mode_labels = load_manual_mode_labels(fpath)
% loads hand annotated labels for unsupervised watershed modes from csv file
% specified by fpath

fID = fopen(fpath);
mode_labels = textscan(fID,'%s %s %s','Delimiter',',');
mode_labels = cat(2,mode_labels{:});
mode_labels(1,:)=[];
mode_labels = cellfun(@(ml) sprintf('%s - %s [%s]',ml{2},ml{3},ml{1}),...
    num2cell(mode_labels,2), 'UniformOutput', false);

