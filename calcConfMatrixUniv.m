function [TNtot, TPtot, FPtot, FNtot, FDtot] = calcConfMatrixUniv(p, corrVector, signThreshold, corrThresh)

Pf = (p < signThreshold);

TP=0;
TN=0;
FP=0;
FN=0;

nVars = size(p, 2);

for i=1:nVars
    if(abs(corrVector(i)) < corrThresh && Pf(i) == 0)
        TN=TN+1;
    elseif(abs(corrVector(i)) > corrThresh && Pf(i) == 1)
        TP=TP+1;
    elseif(abs(corrVector(i)) > corrThresh && Pf(i) == 0)
        FN=FN+1;
    elseif(abs(corrVector(i)) < corrThresh && Pf(i) == 1)
        FP=FP+1;
   end
end

TNtot = TN/(FP+TN);
TPtot = TP/(TP+FN);
FPtot = FP/(FP+TN);
FNtot = FN/(TP+FN);
FDtot = FP/(TP+FP);

end