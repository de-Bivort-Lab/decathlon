function [sigma,z_thresh] = fit_tsne_z_logspeed_gmm(z_speed,plot_bool)
% fit two-component gaussian mixture model to t-SNE trajectory log speed data,
% find speed threshold for parsing data, and (optionally) plot results

% convert speed to log speed and sub-sample the data
log_z = log(cat(1,z_speed{:}));
log_z(isnan(log_z)|isinf(log_z))=[];
sample_sz = 50000;
sample_sz(sample_sz>numel(log_z)) = numel(log_z);
log_z_sample = log_z(randperm(numel(log_z),sample_sz));

% fit two-component gmm to log speed
opts = statset;
opts.MaxIter = 3000;
mdl = fitgmdist(log_z_sample,2,'Options',opts);

% estimate gaussian pdfs based on model
x = linspace(-2,10,200);
[~,centera]=hist(log_z,x);
gaussPdf= pdf(mdl,centera');
gaussPdfi = NaN(numel(x),mdl.NumComponents);
sigma = NaN(mdl.NumComponents,1);
mu = NaN(mdl.NumComponents,1);
weight = NaN(mdl.NumComponents,1);
for n = 1:mdl.NumComponents
  mu(n)          = mdl.mu(n);
  sigma(n)       = sqrt(mdl.Sigma(1,1,n));
  weight(n)      = mdl.ComponentProportion(n);
  gaussPdfi(:,n) = weight(n)*normpdf(centera',mu(n),sigma(n));
end

% sort modes such that low speed mode comes first
[~,p] = sort(mu);
mu = mu(p);
sigma = sigma(p);
weight = weight(p);
gaussPdfi = gaussPdfi(:,p);

% find threshold
[~,mu1_idx] = min(abs(centera-mu(1)));
[~,mu2_idx] = min(abs(centera-mu(2)));
[~,intersect_idx] = min(abs(diff(gaussPdfi(mu1_idx:mu2_idx,:),1,2)));
intersect_idx = intersect_idx + mu1_idx -1;
log_z_thresh = centera(intersect_idx);
z_thresh = exp(log_z_thresh);

% generate optional plot
if plot_bool
    figure;
    hold on;
    hh=histogram(log_z,x,'Normalization','pdf','FaceColor',[0.1 0.1 .3],'FaceAlpha',1,'EdgeColor','none');
    hpdf=plot(centera', gaussPdf, 'Color',[.5 .5 1], 'linewidth', 2);
    hg1=plot(centera', gaussPdfi(:,1),'-m','linewidth', 2);
    hg2=plot(centera', gaussPdfi(:,2),'Color',[0 .75 0],'linewidth', 2);
    ylim = get(gca,'YLim');
    plot([log_z_thresh log_z_thresh],[0 ylim(2)],'r--');
    hold off;
    legend([hpdf;hg1;hg2;hh],{'pdf';'component 1';'component 2';'histogram'});
    title_str = {'embedded log-speed GMM';...
        sprintf('(\\mu_{1}=%0.1f, \\sigma_{1}=%0.1f;  \\mu_{2}=%0.1f, \\sigma_{2}=%0.1f)',...
        mu(1),sigma(1),mu(2),sigma(2))};
    title(title_str);
    ylabel('PDF');
    xlabel('log(speed)');
end