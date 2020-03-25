function [coeffs,r,p] = bootstrap_linmodel(x,y)

coeffs = NaN(size(x,1),size(x,2)+1);
y_true = NaN(size(y));
y_pred = NaN(size(y));

% compute observed p-value
for i=1:size(x,1)
    % leave one out of data
    f = true(size(y));
    f(i) = false;
    mdl = fitlm(x(f,:),y(f));
    
    % define test set
    y_true(i) = y(~f);
    y_pred(i) = predict(mdl,x(~f,:));
end
[r,p] = corr(y_true,y_pred);

