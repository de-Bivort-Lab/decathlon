function [thresh, class_means, sumd] = dec_kthresh_distribution(data,varargin)

% exclude NaN and inf values
data(data<log(1E-4)) = [];
data(isnan(data)) = [];
data(isinf(data)) = [];
nreps = 10;
k=2;
for i=1:numel(varargin)
    switch i
        case 1, nreps = varargin{i};
        case 2, k = varargin{i};
    end
end
if ~isreal(data)
    data = real(data);
end

thresh = NaN(nreps,1);
sumd = NaN(nreps,k);
for i=1:nreps
    if numel(data) > 100000
        sample = data(randperm(numel(data),100000));
    else
        sample = data;
    end

    % perform k-means clustering
    if numel(sample)>2
        [idx,class_means,sumd(i,:)] = kmeans(sample,k,'MaxIter',1000);
        idx = [min(sample(idx==1)) max(sample(idx==1)) min(sample(idx==2)) max(sample(idx==2))];
        class_boundaries = sort(idx);
        class_means = sort(class_means);
        thresh(i) = mean(class_boundaries(2:3));
    else
        thresh(i) = [];
        class_means = [];
    end
end

sumd = nanstd(thresh);
thresh = nanmedian(thresh);
if thresh < 2
%     ah = subplot(2,1,1);
%     hold on;
%     plot(x,y);
%     plot([thresh thresh],[0 max(y)*1.2],'r--');
%     set(ah,'XLim',[-.5 4]);
%     ah = subplot(2,1,2);
%     plot(x,dy2);
%     set(ah,'XLim',[-.5 4]);
%     disp('stop');
end