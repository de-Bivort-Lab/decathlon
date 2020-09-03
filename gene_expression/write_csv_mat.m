function write_csv_mat(fname,data,grp_idx,grp_name,k_ids)

fid = fopen(fname,'W+');
format_spec = ['%s' repmat(',%i',1,size(data,2)) '\n'];
col_names = cellfun(@(s,n) repmat({[', ' s]},1,n), ...
    grp_name, num2cell(cellfun(@numel,grp_idx)), 'UniformOutput', false);
col_names = cat(2,col_names{:});
header = ['kegg_gid' cat(2,col_names{:}), '\n'];
fprintf(fid,header);
cellfun(@(s,p) fprintf(fid,format_spec,s,p), k_ids, num2cell(data,2));
fclose(fid);