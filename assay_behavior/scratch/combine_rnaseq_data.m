% add D1 sequencing data to D2/D3 RNAseq struct
% r1 = D1 RNAseq
% r2 = D2 RNAseq

% add extra index to add D1 data to D2/D3
D_seq = repmat(struct('data',[],'geneID',[],'ID',[]),3,1);
D_seq(1).data = r1.reads';
D_seq(1).ID = (1:192)';
D_seq(1).geneID = r1.geneID;

% D2
reads = r2(1).raw_reads;
data = NaN(size(reads,1),192);
data(:,r2(1).fly_ids-192) = reads;
D_seq(2).data = data';
D_seq(2).geneID = r2(1).geneIDs;
D_seq(2).ID = (193:384)';

% D3 (WAITING ON IDs FROM JULIEN)
% reads = r2(1).raw_reads;
% data = NaN(size(reads,1),192);
% data(:,r2(1).fly_ids-192) = reads;
% D_seq(3).data = data';
% D_seq(3).geneID = r2(1).geneIDs;
% D_seq(3).ID = 193:384;