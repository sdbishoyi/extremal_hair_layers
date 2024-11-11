function [scalar] = SphericalHarmonics_HighEll(ell,m,theta,phi)

euler = 2.7182818284590452353602874713527;

% Computes scalar harmonics to high ell, so needs recursion relation
% Code is vectorized so that theta and phi can be vectors. 
% this code based on Numerical Recipies 

%first compute P_lm where x=cos(theta)

x = cos(theta);
m_sign=1;
if(m>ell)
    scalar=0;
    return
end
if(m<0)
    m_sign=-1;
    m=abs(m);
end
pmm = 1;
if(m>0)
    somx2 = sqrt((1-x).*(1+x));
    fact  = 1.0;
    for i =1:m
       pmm  = -pmm*fact*somx2;
       fact =  fact+2;
    end
end
if(ell==m)
   plgndr = pmm;
else
   pmmp1 = x.*(2*m+1).*pmm;
   if(ell == (m+1))
        plgndr = pmmp1;
   else
        for ll = m+2:ell
           pll = (x.*(2*ll-1).*pmmp1-(ll+m-1).*pmm)./(ll-m);
           pmm   = pmmp1;
           pmmp1 = pll;
        end
        plgndr = pll;
   end
end

%Some factorials needed for norming constant. ifactor is (l+m)!/(l-m)!.

ifactor=1.0;
for j=ell-m+1:ell+m
    ifactor = j*ifactor;
end
const  = 0.5*(2*ell+1)/ifactor;
const  = const/(2*pi);
const  = sqrt(const);
scalar=const.*plgndr.*euler.^(sqrt(-1).*m.*phi);

%relation for m<0
if(m_sign<0)
    scalar = (-1)^m*conj(scalar);
end
return


