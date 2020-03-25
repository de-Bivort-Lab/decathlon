%%

fDir = autoDir;
if ~iscell(fDir)
   fDir = {fDir}; 
end
fPaths = cellfun(@(d) getHiddenMatDir(d), fDir, 'UniformOutput',false);
fPaths = cat(2,fPaths{:});
fPaths = fPaths';
options = {'Raw',true,'Regress',false,...
    'Slide',false,'Bootstrap',false,'Save',true,'Zip',false};


%%

for i=1:numel(fPaths)
    fprintf('processing file %i of %i\n',i,numel(fPaths));
    [tmp_dir,~,~] = fileparts(fPaths{i});
    options = {'Dir',[tmp_dir '\'],'Raw',true,'Regress',false,...
        'Slide',true,'Bootstrap',false,'Save',true,'Zip',false};
    analyze_multiFile(options{:});
end
