function cmap = logjet_cmap

cmap_base = jet;
idx = 1:7:size(cmap_base,1);
log_cmap = fliplr(log(cmap_base(idx,:)+1)')';
cmap = exp(interp1(log(idx),log_cmap,linspace(0,max(log(idx)),256)))-1;
cmap = fliplr(cmap')';