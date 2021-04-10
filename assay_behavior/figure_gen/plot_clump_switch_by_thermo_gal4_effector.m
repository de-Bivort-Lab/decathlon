% Create scatter plots and correlation bar plots for Y-Maze turn clumpiness and 
% switchiness from the thermogenetic gal4 LDM screen. 
% Data input: D_thermo (from decathlon_behavior_data.mat)

% Define the INPUT DIRECTORY
fdir = uigetdir(pwd,'Select decathlon_paper_data directory');
D_thermo = load_decathlon_structs(fdir,'D_thermo');

%%
figure;
colors = {[177 128 255]./255;[102 102 102]./255;[84 147 96]./255};
nplots = 6;
ax = subplot_array(nplots);
delete(ax(3));


% filter out most extreme points
all_effs = D_thermo.effector;
unique_effs = unique(all_effs);
unique_effs = unique_effs([2 1 3]);

x_23 =nanzscore(cellfun(@(d)nanmean(d(1,:)),D_thermo.clumpiness.data));
y_23 =nanzscore(cellfun(@(d)nanmean(d(1,:)),D_thermo.switchiness.data));
x_33 =nanzscore(cellfun(@(d)nanmean(d(3,:)),D_thermo.clumpiness.data));
y_33 =nanzscore(cellfun(@(d)nanmean(d(3,:)),D_thermo.switchiness.data));

% calculate projection of points onto the vector x=1, y=1
d_pts = NaN(numel(x_23),2);
d_pts(:,1) = ([x_23 y_23]*[1;1])./([1 1]*[1;1]);
d_pts(:,2) = ([x_33 y_33]*[1;1])./([1 1]*[1;1]);
d_eff = NaN(numel(unique_effs),2);
r_eff = NaN(numel(unique_effs),2);
ci_d = NaN(3,6);
ci_r = NaN(3,6);
for i=1:numel(unique_effs)
    % group data by effector type
    eff_idx = strcmpi(all_effs,unique_effs{i});
    
    % plot data at low temp
    hold(ax(1),'on');
    pretty_scatter(x_23(eff_idx),y_23(eff_idx),colors{i},'Parent',ax(1));
    xlabel(ax(1),'Clumpiness');
    ylabel(ax(1),'Switchiness');
    r = corr([x_23 y_23],'Type','Spearman','rows','pairwise');
    title(ax(1),sprintf('Temp = 23C; r = %0.2f',r(1,2)));
    dd = abs(d_pts(eff_idx,1));
    d_eff(i,1) = mean(dd);
    ci_d(1:2,i) = bootstrap_mean_CI(dd,0.05,1000);
    
    tmp = [x_23(eff_idx) y_23(eff_idx)];
    tmp(any(isnan(tmp),2),:) = [];
    bs_r = bootstrap_rvalues(tmp,1000);
    r_eff(i,1) = mean(bs_r(1,2,:));
    ci_r(1:2,i) = bootstrap_mean_CI(squeeze(bs_r(1,2,:)),0.05,1000);
    
    % plot data at high temp
    hold(ax(2),'on');
    pretty_scatter(x_33(eff_idx),y_33(eff_idx),colors{i},'Parent',ax(2));
    xlabel(ax(2),'Clumpiness');
    ylabel(ax(2),'Switchiness');
    r = corr([x_33 y_33],'Type','Spearman','rows','pairwise');
    title(ax(2),sprintf('Temp = 33C; r = %0.2f',r(1,2)));
    dd = abs(d_pts(eff_idx,2));
    d_eff(i,2) = mean(dd);
    ci_d(1:2,i+3) = bootstrap_mean_CI(dd,0.05,1000);
    
    tmp = [x_33(eff_idx) y_33(eff_idx)];
    tmp(any(isnan(tmp),2),:) = [];
    bs_r = bootstrap_rvalues(tmp,1000);
    r_eff(i,2) = mean(bs_r(1,2,:));
    ci_r(1:2,i+3) = bootstrap_mean_CI(squeeze(bs_r(1,2,:)),0.05,1000);
end

axis(ax(1:2),'equal');
rr = 3;
plot([-rr rr],[-rr rr],'k--','Parent',ax(1));
plot([-rr rr],[-rr rr],'k--','Parent',ax(2));
for i=1:2
   set(ax(i),'XLim',[-rr rr],'YLim',[-rr rr]); 
end

% plot r-values
bh = bar(r_eff');
ngroups = size(r_eff, 2);
nbars = size(r_eff, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
x = NaN(3,2);
for i=1:numel(bh)
    bh(i).FaceColor = colors{i};
    x(i,1:2) = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
end
x = [x(:)'; x(:)'; NaN(1,numel(x))];
hold on;
plot(x(:), ci_r(:), 'k-','LineWidth',1);
set(ax(6),'YLim',[-1 1],'XTickLabel',{'23C';'29C'});
ylabel('correlation');
legend(bh,{'Shi';'Iso';'Trp'},'Location','SouthEast');

hold(ax(4),'on');
PCARegressionCI([x_23,y_23],ax(4));
xlabel(ax(4),'Clumpiness');
ylabel(ax(4),'Switchiness');

hold(ax(5),'on');
PCARegressionCI([x_33,y_33],ax(5));
xlabel(ax(5),'Clumpiness');
ylabel(ax(5),'Switchiness');

% % generate confidence intervals for temp 23
% hold(ax(4),'on');
% xpts = linspace(-rr,rr,1000);
% nreps = 1000;
% ypts = zeros(nreps,numel(xpts));
% ci95 = zeros(2,numel(xpts));
% for i=1:nreps
%     idx = randi(numel(x_23),[numel(x_23) 1]);
%     b = polyfit(x_23(idx),y_23(idx),1);
%     ypts(i,:) = xpts.*b(1) + b(2);
% end
% ci95(1,:) = prctile(ypts,97.5);
% ci95(2,:) = prctile(ypts,2.5);
% vx = [xpts fliplr(xpts)];
% vy = [ci95(2,:) fliplr(ci95(1,:))];
% patch('XData',vx(:),'YData',vy(:),'FaceColor',[.75 .75 .75],...
%     'EdgeColor','none','Parent',ax(4));
% b = polyfit(x_23,y_23,1);
% plot(xpts,xpts.*b(1) + b(2),'k--','Parent',ax(4));
% pretty_scatter(x_23,y_23,'k','Parent',ax(4));
% axis('equal');
% set(ax(4),'XLim',[-rr rr],'YLim',[-rr rr]);
% xlabel(ax(4),'Clumpiness');
% ylabel(ax(4),'Switchiness');
% 
% % generate confidence intervals for temp 33
% mask = isnan(x_33) | isnan(y_33);
% x_33(mask) = [];
% y_33(mask) = [];
% hold(ax(5),'on');
% xpts = linspace(-rr,rr,1000);
% nreps = 1000;
% ypts = zeros(nreps,numel(xpts));
% ci95 = zeros(2,numel(xpts));
% for i=1:nreps
%     idx = randi(numel(x_33),[numel(x_33) 1]);
%     b = polyfit(x_33(idx),y_33(idx),1);
%     ypts(i,:) = xpts.*b(1) + b(2);
% end
% ci95(1,:) = prctile(ypts,97.5);
% ci95(2,:) = prctile(ypts,2.5);
% vx = [xpts fliplr(xpts)];
% vy = [ci95(2,:) fliplr(ci95(1,:))];
% patch('XData',vx(:),'YData',vy(:),'FaceColor',[.75 .75 .75],...
%     'EdgeColor','none','Parent',ax(5));
% b = polyfit(x_33,y_33,1);
% plot(xpts,xpts.*b(1) + b(2),'k--','Parent',ax(5));
% pretty_scatter(x_33,y_33,'k','Parent',ax(5));
% axis('equal');
% set(ax(5),'XLim',[-rr rr],'YLim',[-rr rr]);
% xlabel(ax(5),'Clumpiness');
% ylabel(ax(5),'Switchiness');

ah = findall(gcf,'Type','axes');
for i=1:numel(ah)
   ah(i).Units = 'inches';
   ah(i).Position(3:4) = 0.8229;
end