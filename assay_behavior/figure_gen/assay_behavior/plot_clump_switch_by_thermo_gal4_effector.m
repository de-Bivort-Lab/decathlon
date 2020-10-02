% Create scatter plots and correlation bar plots for Y-Maze turn clumpiness and 
% switchiness from the thermogenetic gal4 LDM screen. 
% Data input: D_thermo (from decathlon_behavior_data.mat)

% Define the INPUT DIRECTORY
fdir = uigetdir(pwd,'Select decathlon_paper_data directory');
D_thermo = load_decathlon_structs(fdir,'D_thermo');

figure;
colors = {[177 128 255]./255;[102 102 102]./255;[84 147 96]./255};


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
    pretty_scatter(x_23(eff_idx),y_23(eff_idx),colors{i},'MarkerSize',3,'Parent',ax(1));
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
    pretty_scatter(x_33(eff_idx),y_33(eff_idx),colors{i},'MarkerSize',3,'Parent',ax(2));
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
ax(4) = subplot(2,nplots,3);
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
set(ax(4),'YLim',[-1 1],'XTickLabel',{'23C';'29C'});
ylabel('correlation');
legend(bh,{'Shi';'Iso';'Trp'},'Location','SouthEast');

% generate confidence intervals
tmp = [x_23 y_23];
tmp(any(isnan(tmp),2),:) = [];
PCARegressionCI([x_23,y_23],subplot(2,nplots,4),'XLim',[-rr rr],'YLim',[-rr rr]);
axis('equal');
set(gca,'XLim',[-rr rr],'YLim',[-rr rr]);

tmp = [x_33 y_33];
tmp(any(isnan(tmp),2),:) = [];
PCARegressionCI(tmp,subplot(2,nplots,5),'XLim',[-rr rr],'YLim',[-rr rr]);
axis('equal');
set(gca,'XLim',[-rr rr],'YLim',[-rr rr]);