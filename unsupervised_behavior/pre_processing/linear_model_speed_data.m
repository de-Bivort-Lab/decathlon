%% use als to fill data

s = cat(1,all_speed{1:63});
w = [cat(1,all_wavelets{1:63}) s];

%%

witheld_pts = rand([size(w,1) 1])>0.75;
w_miss = w(witheld_pts,end);
w(witheld_pts) = NaN;

opts = statset;
opts.Display = 'iter';
[coefs,score,~,~,~,mu] = pca(w,'algorithm','als','Options',opts);
w_fill = score*coefs' + mu;

%% polyfit

w = zscore(w);
mdl = fitlm(w(~witheld_pts,1:end-1),w(~witheld_pts,end),'linear');

x = w(witheld_pts,1:end-1);
y = mdl.predict(x);
true_y = w(witheld_pts,end);
figure;
pretty_scatter(true_y,y,'k','MarkerSize',1);
figure;
histogram(y-true_y);
resid = y-true_y;
mean(resid.^2)