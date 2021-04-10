function [null_exp, obs_exp, nkeep, plot_handles] = plot_pca_bootstrap(data,nreps,varargin)
% Plot variance explained for PCA performed on bootstrapped observed and 
% shuffled data (null model)

ci=95;
mode='noncummulative';
color = [.5 .5 .5];
max_n = size(data,1);
for i=1:numel(varargin)
   switch i
       case 1, ci = varargin{i};
       case 2, mode = varargin{i};
       case 3, color = varargin{i};
       case 4, max_n = varargin{i};
   end
end

color_null = color - 0.5;
color_null(color_null<0) = 0;

% bootstrap PCA null model by shuffling data matrix and performing PCA
[~, ~, null_exp] = bootstrap_pca_nullmodel(data,nreps,max_n);
[null_ph, null_lh, mu_null_exp, null_ci95] = plot_model(null_exp,color_null,ci,mode);
null_lh.LineStyle = '-';

% bootstrap PCA null model by shuffling data matrix and performing PCA
[~, ~, obs_exp] = bootstrap_pca_observed(data,nreps,max_n);
[obs_ph, obs_lh, mu_obs_exp] = plot_model(obs_exp,color,ci,mode);

plot_handles = {obs_lh;null_lh;obs_ph;null_ph};


nkeep = find(mu_obs_exp>null_ci95(:,1)&mu_obs_exp>1/size(data,2)/10,1,'Last');
if ~isempty(nkeep)
    set(gca,'Ylim',[0.01 100],'Xlim',[1 size(data,2)]);
    plot([nkeep nkeep],[0.01 100],'k--');
    title(sprintf('k = %i',nkeep));
else
    nkeep = 1;
end

% configure plot labels
leg_labels = {'observed mean';'shuffled mean'};
switch mode
    case 'cummulative'
        ylabel({'cummulative';'variance explained'});
        xlabel('num PCs included');
         legend([obs_lh,null_lh],leg_labels,'Location','SouthEast');
    otherwise
        ylabel('variance explained');
        xlabel('PC number');
         legend([obs_lh,null_lh],leg_labels,'Location','NorthEast');
end



function [ph, lh, var_exp, ci95] = plot_model(var_exp,color,ci,mode)

% compute var explained mean and confidence interval
var_exp = cat(2,var_exp{:});
switch mode
    case 'cummulative'
        var_exp = cumsum(var_exp);
        ylim = 100;
    otherwise
        ylim = ceil(max(var_exp(:))*1.1);
end
        
ci95 = NaN(size(var_exp,1),2);
ci95(:,1) = prctile(var_exp,(100-ci)/2,2);
ci95(:,2) = prctile(var_exp,100-(100-ci)/2,2);
var_exp = nanmean(var_exp,2);

% plot mean and CI
npcs = numel(var_exp);
x = 1:npcs;

lh = semilogy(1:npcs,var_exp,'Color',color,'LineStyle','-','LineWidth',1);
hold on;
plot(1:npcs,var_exp,'-','Color',color,'LineWidth',1);
%lh = plot(1:npcs,var_exp,'Color',color,'LineStyle','-','LineWidth',1);
vx = [1 x fliplr(x)];
vy = [ci95(1,2) ci95(:,1)' fliplr(ci95(:,2)')];
ph = patch('XData',vx(:),'YData',vy(:),'FaceColor',color,'FaceAlpha',0.2,...
    'EdgeColor','none','LineWidth',0.5);

%set(gca,'XLim',[1 npcs],'YLim',[0 ylim],'XTick',x,'XTickLabels',1:numel(x));