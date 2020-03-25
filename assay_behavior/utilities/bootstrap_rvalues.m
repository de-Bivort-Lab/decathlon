function corrmats = bootstrap_rvalues(data,nreps,varargin)

max_i = size(data,1);
corrmats = NaN(size(data,2),size(data,2),nreps);
for i=1:nreps
    samples = randi(max_i,[max_i 1]);
    corrmats(:,:,i) = corr(data(samples,:),'Type','Spearman',varargin{:});
end