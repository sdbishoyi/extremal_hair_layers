%%%% self written Newton root solving algorithm using Scott's Sch code as base
%%%% code works but has some discontinuities near horizon and inside horizon
%%%% Gaurav suggested to go from right side and then use the previous r* value
%%%% as guess for next which is used in _Newton.m file
%%%% comment added on 24/01/22


function [rkerr] = RstarToRwithSpin(rs,m,a)
    rkerr = zeros(size(rs));
    rm = m - sqrt(m^2 - a^2) ; 
    rp = m + sqrt(m^2 - a^2) ;
    LENGTHrs = size(rs,1)*size(rs,2);
    
    for k=1:LENGTHrs
        rstar=rs(k);
        if rstar < -1.0e3   %less than -1000, very close to horizon (when m=1)
              error = 'Tortoise coordinate too small'
              pause
        end
        if rstar > 4.0*m   %not near horizon
          r  = rstar;
        
        else                %near horizon
%           rschm2m = 2*m*exp(rstar/2*m - 1);
          r    = 2*rp + 1e-14 ;
        end

        
        for i=1:20   %Newtons method
             if abs(r) < 1.0e-25  %approximation is good enough
                 break
             end
                  funct     = r + (2*m*rp/(rp-rm))*log(abs(r-rp)/2*m) ...
                                - (2*m*rm/(rp-rm))*log(abs(r-rm)/2*m) ...
                                - rstar;
                  fprime    = 1 + 2*m*rp/((rp-rm)*(r-rp)) ...
                                - 2*m*rm/((rp-rm)*(r-rm));
                  r = r - funct/fprime;
%                 r = r - (r + (2*m*rp/(rp-rm))*log((r-rp)/2*m) - (2*m*rm/(rp-rm))*log((r-rm)/2*m) - rstar)/...
%                     (1 + 2*m*rp/((rp-rm)*(r-rp)) - 2*m*rm/((rp-rm)*(r-rm)));
%              if r-rm<0
%                  r=0.5*r;
%              end  
             if r-rp<0
                 r=0.5*r;
             end 
        end
        
        if (r-rp)<1e-15
            r=rp;
        end    
        rkerr(k)=r;
    end
    
    
  
