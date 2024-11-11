function [dlogAbsF_dlogt] = ComputeTailDecayRate(T,F)
%
logT   =log10(T);
logAbsF=log10(abs(F));
%
DlogAbsF = logAbsF(2:end) - logAbsF(1:end-1);
DlogT    = logT(2:end)    - logT(1:end-1);
% analytically the log time spacings should be this when T is equally space
DtTest = transpose((1+(1:length(T)-1))./ (1:length(T)-1));
logDtTest=log10(DtTest);
max(abs(DlogT-logDtTest));
%DlogAbsF = reshape(DlogAbsF,length(DlogAbsF),1);
%DlogT    = reshape(DlogT,length(DlogT),1);
dlogAbsF_dlogt=DlogAbsF./DlogT;