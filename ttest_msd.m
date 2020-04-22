function [ tprob ] = ttest_msd( Amean, Asd, Bmean, Bsd, N)
% TTEST_MSD performs a paired ttest when only the means and standard deviations are available 
% output: tprob = [2-tailed p-value, 1-tailed p-value]
% based on - 
% https://www.mathworks.com/matlabcentral/answers/301296-how-to-do-paired-t-test-with-mean-and-sd
v = 2*N-2;
tval = (Amean-Bmean) / sqrt((Asd^2+Bsd^2)/N); %calcualte T-statistic
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5)); %2-tailed t-distribution
tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;

tprob = 1-[tdist2T(tval,v) tdist1T(tval,v)];
end

