function [output] = PCalc_Continuous(data,effectSizes, sampSizes, signThreshold, nSimSamp, nRepeats)

% If sample size bigger than number of simulated samples adjust it
if max(sampSizes) >= nSimSamp
    disp('Number of simulated samples smaller than maximum of samplesizes to check - increased');
    nSimSamp = max(sampSizes) + 500;
end

% Number of variables
numVars = size(data,2);

% Number of sample and effect sizes
nSampSizes = size(sampSizes,2);
nEffSizes = size(effectSizes, 2);

%Simulation of a new data set based on multivariate normal distribution%

[Samples, correlationMat] = simulateLogNormal(data, 'Estimate', nSimSamp);

% Initialize the output structures
output = cell(1, numVars);

fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,numVars) '\n\n']);

parfor currVar = 1:numVars
    
    storeVar = cell(4, 10);
    
    uncStruct = struct('TP', zeros(nEffSizes, nSampSizes) ,'FP',  zeros(nEffSizes, nSampSizes) ,'TN',  zeros(nEffSizes, nSampSizes), ...
        'FN',  zeros(nEffSizes, nSampSizes), 'FD',  zeros(nEffSizes, nSampSizes), 'STP', zeros(nEffSizes,nSampSizes), 'SFP', zeros(nEffSizes,nSampSizes), ...
        'STN', zeros(nEffSizes, nSampSizes), 'SFN', zeros(nEffSizes, nSampSizes), 'SFD', zeros(nEffSizes, nSampSizes));
    
    bonfStruct = struct('TP', zeros(nEffSizes, nSampSizes) ,'FP',  zeros(nEffSizes, nSampSizes) ,'TN',  zeros(nEffSizes, nSampSizes), ...
        'FN',  zeros(nEffSizes, nSampSizes), 'FD',  zeros(nEffSizes, nSampSizes), 'STP', zeros(nEffSizes,nSampSizes), 'SFP', zeros(nEffSizes,nSampSizes), ...
        'STN', zeros(nEffSizes, nSampSizes), 'SFN', zeros(nEffSizes, nSampSizes), 'SFD', zeros(nEffSizes, nSampSizes));
    bhStruct = struct('TP', zeros(nEffSizes, nSampSizes) ,'FP',  zeros(nEffSizes, nSampSizes) ,'TN',  zeros(nEffSizes, nSampSizes), ...
        'FN',  zeros(nEffSizes, nSampSizes), 'FD',  zeros(nEffSizes, nSampSizes), 'STP', zeros(nEffSizes,nSampSizes), 'SFP', zeros(nEffSizes,nSampSizes), ...
        'STN', zeros(nEffSizes, nSampSizes), 'SFN', zeros(nEffSizes, nSampSizes), 'SFD', zeros(nEffSizes, nSampSizes));
   
    byStruct = struct('TP', zeros(nEffSizes, nSampSizes) ,'FP',  zeros(nEffSizes, nSampSizes) ,'TN',  zeros(nEffSizes, nSampSizes), ...
        'FN',  zeros(nEffSizes, nSampSizes), 'FD',  zeros(nEffSizes, nSampSizes), 'STP', zeros(nEffSizes,nSampSizes), 'SFP', zeros(nEffSizes,nSampSizes), ...
        'STN', zeros(nEffSizes, nSampSizes), 'SFN', zeros(nEffSizes, nSampSizes), 'SFD', zeros(nEffSizes, nSampSizes));
    
    for currEff=1:nEffSizes
        
        b1 = zeros(numVars,1);
        
        b1(currVar) = effectSizes(currEff);
        
        for currSampSize = 1:nSampSizes
            
            % initialize arrays to store the repeats
            multiplerepeats = struct();
            
            multiplerepeats.Results = struct('TP', zeros(1,nRepeats) ,'FP',  zeros(1, nRepeats) ,'TN',  zeros(1, nRepeats), ...
                'FN',  zeros(1, nRepeats), 'FD',  zeros(1, nRepeats));
            
            multiplerepeats.Bonferroni = struct('TP', zeros(1, nRepeats) ,'FP',  zeros(1, nRepeats) ,'TN',  zeros(1, nRepeats), ...
                'FN',  zeros(1, nRepeats), 'FD',  zeros(1, nRepeats));
            
            multiplerepeats.BenjHoch = struct('TP', zeros(1, nRepeats) ,'FP',  zeros(1, nRepeats) ,'TN',  zeros(1, nRepeats), ...
                'FN',  zeros(1, nRepeats), 'FD',  zeros(1, nRepeats));
            
            multiplerepeats.BenjYek = struct('TP', zeros(1, nRepeats) ,'FP',  zeros(1, nRepeats) ,'TN',  zeros(1, nRepeats), ...
                'FN',  zeros(1, nRepeats), 'FD',  zeros(1, nRepeats));
            
            for currRepeat=1:nRepeats
                
                selectIndex = randperm(sampSizes(currSampSize));
                
                SelSamples = Samples(selectIndex, :);

                % UVScaling the data - vectorize with bsxfun
                stDev = std(SelSamples, 1);
                SelSamples = bsxfun(@minus, SelSamples, mean(SelSamples, 1));
                SelSamples = bsxfun(@rdivide, SelSamples, stDev);
                
                noiseLevel = 1;
                noise = noiseLevel*randn(sampSizes(currSampSize),1);
                
                Y = SelSamples(:, currVar)*b1(currVar);
                Y = Y + noise;
                %Y = Uvscaled(Y);
                
                
                p = zeros(1, numVars);
             
                % fitlm if needed
                %for tocheck = 1:numVars
                    %mdl = fitlm(SelSamples(:, tocheck), Y);
                    %betafound =  mdl.Coefficients.Estimate('x1');
                    %betapval = mdl.Coefficients.pValue('x1');
                    %p(tocheck) = betapval;
                %end
                
		       % Using regress for multivariate regression
		       % test
                for i = 1:numVars
                    B = [ones(size(Y,1),1), SelSamples(:,i)];
                    [b, bint , r , rint, stats] = regress(Y,B);
                    %b2(i) = b(2,1);
                    p(i) = stats(3);
                end
                
                pUnc = p;
                pBonf = numVars*p;
                [h, crit_p, adj_ci_cvrg, pBY] = fdr_bh(p, 0.05, 'dep');
                [h, crit_p, adj_ci_cvrg, pBH] = fdr_bh(p, 0.05, 'pdep');
                
                corrVector = correlationMat(:, currVar);
                
                [uncTNTot, uncTPTot, uncFPTot, uncFNTot, uncFDTot] = calcConfMatrixUniv(pUnc, corrVector, signThreshold, 0.8);
                [bonfTNTot, bonfTPTot, bonfFPTot, bonfFNTot, bonfFDTot] = calcConfMatrixUniv(pBonf, corrVector, signThreshold, 0.8);
                [byTNTot, byTPTot, byFPTot, byFNTot, byFDTot] = calcConfMatrixUniv(pBY, corrVector, signThreshold, 0.8);
                [bhTNTot, bhTPTot, bhFPTot, bhFNTot, bhFDTot] = calcConfMatrixUniv(pBH, corrVector, signThreshold, 0.8);
                
                multiplerepeats.noCorrection.TP(currRepeat) = uncTPTot;
                multiplerepeats.noCorrection.FP(currRepeat) = uncFPTot;
                multiplerepeats.noCorrection.TN(currRepeat) = uncTNTot;
                multiplerepeats.noCorrection.FN(currRepeat) = uncFNTot;
                multiplerepeats.noCorrection.FD(currRepeat) = uncFDTot;
                
                multiplerepeats.Bonferroni.TP(currRepeat) = bonfTPTot;
                multiplerepeats.Bonferroni.FP(currRepeat) = bonfFPTot;
                multiplerepeats.Bonferroni.TN(currRepeat) = bonfTNTot;
                multiplerepeats.Bonferroni.FN(currRepeat) = bonfFNTot;
                multiplerepeats.Bonferroni.FD(currRepeat) = bonfFDTot;
                
                multiplerepeats.BenjHoch.TP(currRepeat) = bhTPTot;
                multiplerepeats.BenjHoch.FP(currRepeat) = bhFPTot;
                multiplerepeats.BenjHoch.TN(currRepeat) = bhTNTot;
                multiplerepeats.BenjHoch.FN(currRepeat) = bhFNTot;
                multiplerepeats.BenjHoch.FD(currRepeat) = bhFDTot;
                
                multiplerepeats.BenjYek.TP(currRepeat) = byTPTot;
                multiplerepeats.BenjYek.FP(currRepeat) = byFPTot;
                multiplerepeats.BenjYek.TN(currRepeat) = byTNTot;
                multiplerepeats.BenjYek.FN(currRepeat) = byFNTot;
                multiplerepeats.BenjYek.FD(currRepeat) = byFDTot;
              
                
            end
                
                stats = fieldnames(multiplerepeats.Bonferroni);
                for nstat = 1:5
                    currstat = stats{nstat};
                    
                    uncStruct.(currstat)(currEff, currSampSize) =  mean(multiplerepeats.noCorrection.(currstat));
                    uncStruct.(strcat('S', currstat))(currEff, currSampSize) =  std(multiplerepeats.noCorrection.(currstat));
                    bonfStruct.(currstat)(currEff, currSampSize) =  mean(multiplerepeats.Bonferroni.(currstat));
                    bonfStruct.(strcat('S', currstat))(currEff, currSampSize) =  std(multiplerepeats.Bonferroni.(currstat));
                    byStruct.(currstat)(currEff, currSampSize) =  mean(multiplerepeats.BenjYek.(currstat));
                    byStruct.(strcat('S', currstat))(currEff, currSampSize) =  std(multiplerepeats.BenjYek.(currstat));
                    bhStruct.(currstat)(currEff, currSampSize) =  mean(multiplerepeats.BenjHoch.(currstat));
                    bhStruct.(strcat('S', currstat))(currEff, currSampSize) =  std(multiplerepeats.BenjHoch.(currstat));
                    
                end     
         end
                       
    end
    
    stats = fieldnames(uncStruct);
    
    for i = 1:numel(stats)
        storeVar{1, i} = uncStruct.(stats{i});
        storeVar{2, i} = bonfStruct.(stats{i});
        storeVar{3, i} = bhStruct.(stats{i});
        storeVar{4, i} = byStruct.(stats{i});
        
    end
    
    output{1, currVar} = storeVar;
    
    fprintf('\b|\n');
end


