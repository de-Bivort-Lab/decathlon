function [individual_modes,unique_modes,mode_pdfs,idxMap] = ...
                                mode_from_embeddingValues(density,xx,z)

LL = watershed(-density,8);
unique_modes = unique(LL);
LL(density<1E-6)=0;
vals = round((cat(1,z{:}) + max(xx))*length(xx)/(2*max(xx)));
vals(vals==0)=1;
vals_idx = sub2ind(size(LL),vals(:,2),vals(:,1));
watershedValues = LL(vals_idx);
mode_cts  = histc(watershedValues,unique_modes);

minTemplateLength = 10;
rm_mode_idx = unique_modes(mode_cts < minTemplateLength);
modes = watershedValues;
idxMap = LL;
for i=1:length(rm_mode_idx)
    idxMap(LL == rm_mode_idx(i)) = 0;
    idxMap(LL>rm_mode_idx(i)) = idxMap(LL>rm_mode_idx(i)) - 1;
    modes(watershedValues == rm_mode_idx(i)) = 0;
    modes(watershedValues>rm_mode_idx(i)) = modes(watershedValues>rm_mode_idx(i)) - 1;
end
unique_modes = 1:(numel(unique_modes)-numel(rm_mode_idx)-1);

z_lengths = cumsum([0;cellfun(@numel,z)./2]);
individual_modes = cell(numel(z),1);
for i=1:numel(z)
    individual_modes{i} = modes(z_lengths(i)+1:z_lengths(i+1));
end
mode_cts = cellfun(@(m) histc(m,unique_modes), individual_modes, 'UniformOutput', false);
mode_pdfs = cellfun(@(m) m./sum(m(:)), mode_cts, 'UniformOutput', false);
mode_pdfs = cat(2,mode_pdfs{:})';