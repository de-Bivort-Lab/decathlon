function out=decathlonpcaresamp(data)

numBSReps=30;
numFlies=size(data,1);
numDims=size(data,2);
maxPCs=25;
maxPCs(maxPCs>numDims) = numDims;

matchParams.visBool=0;              
matchParams.shuffleTries=100000;    
matchParams.maxPCs=maxPCs;             
matchParams.swapCount=4;  

scoreMatrices=zeros(numDims,maxPCs,numBSReps);

[coeffFull,~,latentFull]=pca(data);

parfor i=1:numBSReps
    
    fprintf(1,'\t replicate #%4i out of %4i\n',i,numBSReps);
    which=randi(numFlies,numFlies,1);
    dataTemp=data(which,:);
    [coeffTemp,~,~]=pca(dataTemp);
    
    matched=matchCoeffMatrices(coeffFull,coeffTemp,latentFull,matchParams);
    
    scoreMatrices(:,:,i)=matched.bestCoeffMatrix(:,1:maxPCs);
    
end

out.scoreMatrices=scoreMatrices;
out.zMatrix=mean(scoreMatrices,3)./std(scoreMatrices,0,3);