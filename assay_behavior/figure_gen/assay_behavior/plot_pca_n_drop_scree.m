function plot_pca_n_drop_scree(data,max_reps)

% iteratively increase the number of dropped features
n_dropped = 1:size(data,2)-1;
intersection = NaN(numel(n_dropped),1);
for i=n_dropped
    % print update
    fprintf('dropping %i features up to %i\n',i,n_dropped(end));
    
    % compute for observed data
%     obs_var_exp = NaN(max_reps,size(data,2)-n_dropped(i));
%     for j=1:max_reps
%         obs_var_exp(j,:) = mean(drop_n_features_pca(data,n_dropped(i),5,false));
%     end
    obs_var_exp = drop_n_features_pca(data,n_dropped(i),max_reps,false);
   
    % compute for shuffled data
    shuf_var_exp = drop_n_features_pca(data,n_dropped(i),max_reps,true);
    
    % find intersection
    below_thresh = obs_var_exp >= prctile(shuf_var_exp,2.5);
    tmp_int = cellfun(@(int) find(int,1,'Last'), ...
        num2cell(below_thresh,2),'UniformOutput',false);
    tmp_int(cellfun(@isempty,tmp_int)) = {size(data,2) - n_dropped(i)};
    tmp_int = cat(1,tmp_int{:});
    intersection(i) = mean(tmp_int);
%     tmp_int = find(mean(obs_var_exp) <= prctile(shuf_var_exp,2.5),1);
%     if isempty(tmp_int)
%         tmp_int = size(data,2) - n_dropped(i);
%     end
%     intersection(i) = tmp_int;
end


% plot results
plot(n_dropped,intersection,'k','LineWidth',1);
xlabel('metrics dropped');
ylabel({'effective';'dimensionality'});
set(gca,'XLim',[1 n_dropped(end)],'YLim',[1 n_dropped(end)]);

hold on;
plot([1 n_dropped(end)],[n_dropped(end) 1],'k--');
