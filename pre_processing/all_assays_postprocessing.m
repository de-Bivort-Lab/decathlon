%%

fDir = autoDir;
fPaths = recursiveSearch(fDir,'ext','.mat');

options = struct('disable',0,'handedness',1,'bouts',1,'bootstrap',0,...
    'slide',0,'regress',0,'areathresh',1,'save',1,'raw',{'speed'});

%%

fprintf('\n');
for i=1:numel(fPaths)
   
   fprintf('processing file %i of %i\n',i,numel(fPaths))
   load(fPaths{i});
   if strcmpi(expmt.meta.name,'Circadian') && expmt.parameters.mm_per_pix==1
       expmt.parameters.mm_per_pix = 0.085;
   elseif ~strcmpi(expmt.meta.name,'Olfaction')
       if expmt.parameters.mm_per_pix==1 && expmt.meta.labels_table.Box(1) == 3
           expmt.parameters.mm_per_pix = 0.3132;
       elseif expmt.parameters.mm_per_pix==1 && expmt.meta.labels_table.Box(1) == 1
           expmt.parameters.mm_per_pix = 0.3302;
       end
   end

   expmt.meta.options = options;
   expmt = autoAnalyze(expmt);
end