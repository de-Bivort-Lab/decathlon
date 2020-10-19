function D = preprocess_rnaseq_data(D, min_tot_reads, min_rpm)
% Pre-process decathlon RNAseq data:
% 1. Sort gene order across decathlon batches
% 2. Filter out flies with too few reads
% 3. Normalize to reads per million (RPM)
% 4. Quantile normalize across individual read count distributions

% use gene names to align reads across decathlon batches
D = pair_rnaseq(D);

% filter read data, update individual fly IDs, and normalize to RPM
reads = cell(numel(D),1);
for i=1:numel(reads)
    reads{i} = D(i).data;
    
    % filter out flies with too few reads
    fly_filt = nansum(reads{i},2) > min_tot_reads;
    reads{i} = reads{i}(fly_filt,:);
    D(i).ID = D(i).ID(fly_filt);
    if i==2
        D(i).ID = D(i).ID - 192;
    end
    
    % normalize to reads per million
    tot_reads = sum(reads{i},2);
    scale_factor = 1E6 ./ tot_reads;
    reads{i} = reads{i} .* repmat(scale_factor,1,size(reads{i},2));
end

% combine reads across batches and quantile normalize reads
all_reads = cat(1,reads{:});
gene_filt = nanmean(all_reads) >= min_rpm;
all_reads = all_reads(:,gene_filt);
q_norm = quantile_normalize(all_reads')';

% partition normalized reads back into batches
for i=1:numel(reads)
    n_rows = [0;cellfun(@(x) size(x,1), reads(1:i))];
    D(i).data = q_norm(sum(n_rows(1:i))+1:sum(n_rows),:);
end

% update paired data struct reads and ID no.
for i=1:numel(D)   
    D(i).geneID = D(i).geneID(gene_filt);
end
