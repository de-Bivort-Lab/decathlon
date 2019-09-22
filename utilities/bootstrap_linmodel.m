function [coeffs,rvals,pvals] = bootstrap_linmodel(x,y,nreps)

coeffs = NaN(nreps,size(x,2)+1);
n_max = size(x,1);
n_train = ceil(n_max*0.7);
rvals = NaN(nreps,1);
pvals = NaN(nreps,1);

for i=1:nreps
    
    % select random subset of data
    train_pts = randperm(n_max,n_train);
    mdl = fitlm(x(train_pts,:),y(train_pts));
    coeffs(i,:) = mdl.Coefficients{:,1};
    
    % define test set
    test_pts = find(~ismember(1:n_max,train_pts));
    y_true = y(test_pts);
    y_pred = predict(mdl,x(test_pts,:));
    
    [r,p] = corr([y_true,y_pred]);
    rvals(i) = r(1,2);
    pvals(i) = p(1,2);
end