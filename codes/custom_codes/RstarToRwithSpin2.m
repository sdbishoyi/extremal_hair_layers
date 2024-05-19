%%%% just a temporary test on how the inbuilt fzero function compares against self written
%%%% Newton root solving algorithm----file needs to be ignored for future implementation ---comment added on 24/01/22

function [rkerr]=  RstarToRwithSpin2(rs,m,a)
    rkerr = zeros(size(rs));
    rm = m - sqrt(m^2 - a^2) ; 
    rp = m + sqrt(m^2 - a^2) ;
    LENGTHrs = size(rs,1)*size(rs,2);
    % rs=[-100:1:100];
    % rs=linspace(-100,100,100);
    rkerr=zeros(size(rs));
    for k=1:LENGTHrs
        rstar=rs(k);
        if rstar > 2*rp   %not near horizon
            rguess    = rstar;
    %           rschm2m = rsch - twom;
        else
            rguess= rp+1e-14;
        end
        
%         if abs(rstar-rp)<1e-3                %near horizon
%     %           rschm2m = twom*exp(rstar/twom - 1);
%               rguess    = rp+1e-15;
%         else
%               rguess    = rp+1e-15;
%         end
%         if r-rm<0
%                  r=0.5*r;
%         end  
%         if r-rp<0
%                  r=0.5*r;
%         end

        fun=@(r)(r + m*log(r^2-2*m*r+a^2) ...
                    +((2*m^2-a^2)/(2*sqrt(m^2-a^2)))*log( (r-m-sqrt(m^2-a^2))/((r-m+sqrt(m^2-a^2))) )...
                    - rstar);
        
%         fun=@(r)(r + (2*m*rp/(rp-rm))*log(abs(r-rp)/2*m) - (2*m*rm/(rp-rm))*log(abs(r-rm)/2*m) - rstar);
        rtemp=fzero(fun,rguess);
%         fprintf('%1.4f \t  %1.4f \n', rtemp, rp);

        if (rtemp-rp)<1e-15
            rtemp=rp;
        end    
        rkerr(k)=rtemp;

    end
    
