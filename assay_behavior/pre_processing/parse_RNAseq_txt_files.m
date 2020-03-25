
RNAseq = repmat(struct('reads',[],'geneIDs',[]),1,1);
fdir = 'D:\decathlon_preprint_code_data_figures\decathlon_paper_data\';
for i=1:numel(RNAseq)
    
    fprintf('READING FILE %i of %i\n',i,numel(RNAseq));
    
    % read in data
    p = sprintf('%sreadcountsallsamples_%i.txt',fdir,i);
    fid = fopen(p);
    data = textscan(fid,'%s','Delimiter','\t','EndOfLine','\r\n');
    data = data{1};
    
    fprintf('Parsing data\n');
    
    % parse header and data
    header_end_idx = find(cellfun(@any,regexp(data,'^(Ben)')),1,'Last');
%     colLabels = data(2:header_end_idx);
    geneID_idx = find(cellfun(@any,regexp(data,'^(FBgn)')));
    if mod(numel(data),header_end_idx)
        data(end-mod(numel(data),header_end_idx)+1:end) = [];
    end
    data = reshape(data,header_end_idx,numel(data)/header_end_idx)';
    RNAseq(i).geneIDs = data(2:end,1);
    col_labels = data(1,2:end)';
    data = data(2:end,2:end);
    data = cellfun(@(s) str2double(s), data(:));
    data = reshape(data,numel(data)/(header_end_idx-1),header_end_idx-1);
    
    fprintf('Sorting data\n');
    
    % parse sample labels into row/col/well
    seq_plate = regexp(col_labels,'(?<=Ben)(\w|[0-9])+','match');
    seq_plate = cat(1,seq_plate{:});
    RNAseq(i).seq_plate = cat(1,seq_plate{:});
    row = regexp(col_labels,'(?<=RNA-)(\w)','match');
    row = cat(1,row{:});
    col = regexp(col_labels,'(?<=RNA-\w)[0-9]+','match');
    col = cat(1,col{:});
    sample_id = cellfun(@(pl,r,c) sprintf('%s-%s-%02d',pl,r,str2double(c)),...
        seq_plate, row, col, 'UniformOutput', false);
    
    % sort samples into ascending well order
    [~,pp] = sort(sample_id);
    RNAseq(i).sampleID = sample_id(pp);
    RNAseq(i).row = row(pp);
    RNAseq(i).col = col(pp);
    RNAseq(i).plate = seq_plate(pp);
    RNAseq(i).raw_reads = data(:,pp);
end

%% Plot total read count for each sample in both files

figure;
for i=1:numel(RNAseq)
    subplot(2,1,i);
    total_reads = sum(RNAseq(i).raw_reads);
    bar(total_reads);
    xlabel('Sample');
    ylabel('Total reads (all genes)');
    title(sprintf('file#%i',i));
    ah = gca;
    set(ah,'XLim',[0 numel(total_reads)+0.5],'TickLength',[0 0],...
        'XTick',1:numel(total_reads),'XTickLabel',RNAseq(i).sampleID,...
        'XTickLabelRotation',90);
    x = [repmat(96:96:96*3,2,1)+0.5; NaN(1,3)];
    y = repmat([0;ah.YLim(2);NaN],3,1);
    hold on;
    plot(x(:),y,'k--','LineWidth',0.5);
end

%% visualize well plate position of low read count / empty samples

well_reads = NaN(12,8,3);
plate = [repmat({'D1'},96,1);repmat({'D3'},96,1);repmat({'E1'},96,1)];
row = repmat('ABCDEFGH',12,3);
col = repmat((1:12)',8,3);
plate_ids = arrayfun(@(p,r,c) sprintf('%s-%s-%02d',p{1},r,c),...
        plate(:), row(:), col(:), 'UniformOutput', false);

sample_id = RNAseq(1).sampleID;
total_reads = sum(RNAseq(1).raw_reads);
for i=1:numel(sample_id)
    idx = find(strcmp(plate_ids,sample_id{i}));
    [ir,ic,ip] = ind2sub(size(well_reads),idx);
    well_reads(ir,ic,ip) = total_reads(i);
end
    
figure;
plate_labels = {'D1';'D3';'E1'};
for i=1:3
    subplot(3,1,i);
    imagesc(well_reads(:,:,i)'>1.5E5);
    colormap('gray');
    set(gca,'XTick',1:size(well_reads,1),'YTick',1:size(well_reads,2),...
        'YTickLabel',num2cell('ABCDEFGH'),'TickLength',[0 0]);
    x = [repmat(1:11,2,1)+.5;NaN(1,11)];
    y = [repmat([.5 12.5]',1,11); NaN(1,11)];
    hold on;
    plot(x(:),y(:),'b-','LineWidth',1);
    y = [repmat(1:7,2,1)+.5;NaN(1,7)];
    x = [repmat([.5 12.5]',1,7); NaN(1,7)];
    plot(x(:),y(:),'b-','LineWidth',1);
    title(sprintf('Plate - %s',plate_labels{i}));
end

%% determine if duplicate geneIDs are numerically identical

[dup_genes_a,idx_b] = ismember(RNAseq(1).geneIDs,RNAseq(2).geneIDs);
idx_b = idx_b(idx_b>0);
x = NaN(sum(dup_genes_a)*282,1);
y = x;
ct = 0;

for i=1:numel(RNAseq(1).geneIDs)
    if ~mod(i,100)
        fprintf('%i\n',i);
    end
    gid = RNAseq(1).geneIDs(i);
    idx = strcmp(geneids,gid);
    if any(idx)
       x(ct+1:ct+282) = tc1(i,:);
       y(ct+1:ct+282) = tc2(idx,:);
       ct = ct+282;
    end
end

x = x(~isnan(x));
y = y(~isnan(y));
scatter(x,y);

%% filter out duplicate genes
    
RNAseq = RNAseq(1);
RNAseq(2) = RNAseq(1);
is_bkiso = true(numel(RNAseq(1).row),1);
is_bkiso(96:191) = false;
f = fieldnames(RNAseq);
for j=1:2
    for i=1:numel(f)
        sz = size(RNAseq(j).(f{i}));
        if sz(1)==282
            if mod(j,2)
                RNAseq(j).(f{i}) = RNAseq(j).(f{i})(is_bkiso,:);
            else
                RNAseq(j).(f{i}) = RNAseq(j).(f{i})(~is_bkiso,:);
            end
        elseif sz(2)==282
            if mod(j,2)
                RNAseq(j).(f{i}) = RNAseq(j).(f{i})(:,is_bkiso);
            else
                RNAseq(j).(f{i}) = RNAseq(j).(f{i})(:,~is_bkiso);
            end
        end
    end
end

%% map sequencing plate IDs to freezing plate IDs

% load freeze to seq map
bk_freeze_map = load('D:\decathlon_data_and_analysis\decathlon 2-2018\meta\D2_freeze_block_to_sequencing_map.mat');
bk_freeze_map = bk_freeze_map.bk_seq_map';
bk_freeze_map = bk_freeze_map(:);

% create plate-row-col label for each sample
plate = [repmat({'D1'},96,1);repmat({'E1'},96,1)];
row = repmat('ABCDEFGH',12,2);
col = repmat((1:12)',8,2);
plate_ids = arrayfun(@(p,r,c) sprintf('%s-%s-%02d',p{1},r,c),...
        plate(:), row(:), col(:), 'UniformOutput', false);

% filter out samples with no entry in freeze map
empty_wells = strcmp(bk_freeze_map,'NaN');
bk_freeze_map(empty_wells) = [];
plate_ids(empty_wells) = [];

% filter out RNAseq samples with no corresponding freezing block ID
sample_ids = RNAseq(1).sampleID;
no_freeze_id = ~ismember(sample_ids,plate_ids);
sample_ids(no_freeze_id) = [];
RNAseq(1).raw_reads(:,no_freeze_id) = [];
RNAseq(1).row(no_freeze_id) = [];
RNAseq(1).col(no_freeze_id) = [];
RNAseq(1).plate(no_freeze_id) = [];
RNAseq(1).sampleID(no_freeze_id) = [];


% filter out freeze IDs with no sequencing sampleID
no_rnaseq_id = ~ismember(plate_ids,sample_ids);
plate_ids(no_rnaseq_id) = [];
bk_freeze_map(no_rnaseq_id) = [];

% parse out freezing block-row-col and create unique identifier
fr_block = regexp(bk_freeze_map,'(?<=_)[A-Z](?=_)','match');
fr_block = cat(1,fr_block{:});
fr_row = regexp(bk_freeze_map,'(?<=_)[A-Z](?=[0-9])','match');
fr_row = cat(1,fr_row{:});
fr_col = regexp(bk_freeze_map,'(?<=_[A-Z])[0-9]+','match');
fr_col = cellfun(@str2double,cat(1,fr_col{:}));
fr_id = cellfun(@(b,r,c) sprintf('%s-%s-%02d',b,r,c),...
    fr_block,fr_row,num2cell(fr_col),'UniformOutput',false);

% assign freezing block meta data to struct
RNAseq(1).freeze = ...
    struct('id',fr_id,'plate',fr_block,'row',fr_row,'col',num2cell(fr_col));


%% build map to freezing block IDs

load('D:\decathlon_data_and_analysis\decathlon 2-2018\meta\D2_freezing_blocks.mat');
D2freezingblocks = D2freezingblocks';
D2freezingblocks = D2freezingblocks(:);
plate = [repmat({'A'},96,1);repmat({'B'},96,1);repmat({'C'},96,1)];
row = repmat('ABCDEFGH',12,3);
col = repmat((1:12)',8,3);
plate_ids = arrayfun(@(p,r,c) sprintf('%s-%s-%02d',p{1},r,c),...
        plate(:), row(:), col(:), 'UniformOutput', false);
    
empty_wells = D2freezingblocks < 1;
D2freezingblocks(empty_wells) = [];
plate_ids(empty_wells) = [];

pp = cellfun(@find,num2cell(str_list_contains(fr_id,plate_ids),2));
RNAseq(1).fly_ids = D2freezingblocks(pp);


