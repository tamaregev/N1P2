function [info] = extractInfolme(tmplme,islm)
%EXTRACTINFOLME extracts only the statistical info from tlme of the 
% LinearMixedModel class (output of fitlme function), without all the data
% which makes it lighter in order to be saved many times.

% input: 
%       tmplme output of fitlme function
%       islm true if tmplme is output of fitlm


if ~exist('islm','var')
    islm=false;
end
% Tamar Regev Jul 18 2019
info.Coefficients = tmplme.Coefficients;
info.CoefficientCovariance = single(tmplme.CoefficientCovariance); 
info.DFE = single(tmplme.DFE);

info.Formula = tmplme.Formula;
info.LogLikelihood = single(tmplme.LogLikelihood);
info.ModelCriterion = tmplme.ModelCriterion;
info.MSE = single(tmplme.MSE);
info.NumObservations = single(tmplme.NumObservations); 
info.Rsquared = tmplme.Rsquared;
info.SSE = single(tmplme.SSE);
info.SSR = single(tmplme.SSR);
info.SST = single(tmplme.SST);
info.anova = anova(tmplme);
if ~islm
    info.FitMethod = tmplme.FitMethod;  
    info.randomEffects = single(randomEffects(tmplme));
    info.fixedEffects = single(fixedEffects(tmplme));
else
    
end

end

