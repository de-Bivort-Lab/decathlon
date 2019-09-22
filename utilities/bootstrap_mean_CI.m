function ci = bootstrap_mean_CI(x,alpha,numReps)

mus = NaN(numReps,size(x,2));
for i=1:numReps
    bs_data = x(randi(size(x,1),[size(x,1) 1]),:);
    mus(i,:) = nanmean(bs_data);
end

ci = NaN(size(x,2),2);
ci(:,1) = prctile(mus,100*alpha/2);
ci(:,2) = prctile(mus,100*(1-alpha/2));