function data = read_csv_mat(fpath)

% open file and read raw text
fid = fopen(fpath);
raw = textscan(fid,'%s','delimiter',',');
raw = raw{1};
fclose(fid);

% find number of columns to initialize format spec
ncol = find(~cellfun(@(s) any(strfind(s,'"')), raw),1) - 2;
header = ['%q' repmat(' %q',1,ncol-1)];
fspec = ['%q' repmat(' %s',1,ncol-1)];

% reopen file and read with correct format spec
fid = fopen(fpath);
head_txt = textscan(fid,header,1,'delimiter',',');
mat_txt = textscan(fid,fspec,'delimiter',',');
mat_txt = [head_txt; cat(2,mat_txt{:})];

% assign data to struct
data.rows = mat_txt(2:end,1);
data.cols = cat(1,mat_txt{1,2:end});
data.data = cellfun(@str2double,mat_txt(2:end,2:end));


