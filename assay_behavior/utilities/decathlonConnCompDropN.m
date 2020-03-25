function out=decathlonConnCompDropN(data)

numReps=20;
numDims=size(data,2);

out=zeros(numDims,numDims+1);

h=waitbar(0,'progress');

for i=0:numDims-2
    
    spectrum=[];
    for j=1:numReps
        whichToDrop=randperm(numDims);
        data_temp=data;
        data_temp(:,whichToDrop(1:i))=[];
        spectrum=[spectrum;decathlonConnCompSweep(data_temp,100)];
    end
    
    out(i+1,1:size(spectrum,2))=mean(spectrum);
    waitbar((i+1)/(numDims+1),h);
    
end

close(h);