function [simData, corrMatrix] = simulateLogNormal(data, covType, nSamples)

offset = abs(min(min(data))) + 1; 

offData = data + offset;

logData = log(offData);

meansLog = mean(logData);

if strcmp(covType, 'Estimate') 
    covLog=cov(logData);
    
elseif strcmp(covType, 'Diagonal') 
     covLog = diag(ones(size(data,2),1));
     covLog(logical(eye(size(covLog)))) = var(logData);
else
    error('Unknown Covariance type');
end

simData = mvnrnd(meansLog',covLog,nSamples);

simData = exp(simData);

simData =  simData - offset;

% Set to 0 negative values 
simData(simData < 0) = 0;

corrMatrix = corr(simData);

end

