function cmap = density_cmap

colors = [255 255 255; 149 174 181; 64 64 89; 0 0 0]./255;
cmap = interp1([0 45 128 205],colors,0:205);
