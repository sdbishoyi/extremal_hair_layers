function [Rrsn,Rschm2m] = RstarToR_ERN(x,m)

%Uses newton's method to turn tortoise coordinates to standard Reissner-Nor
%-strom coordinates
%rstar = read in rstar coordinates, rrsn = RN r coordinate..but
%computed iteratively using newtons method
%rschm2m + 2M ~ rsch... so rschm2m near 0 at horizon instead of 2m,
%we actually update rschm2m and add 2m to this number
%rsn = rstar (newton)

      Rrsn = zeros(size(x));
      Rschm2m = zeros(size(x));
      LENGTHx = size(x,1)*size(x,2);
      twom = m ;
      for k = 1:LENGTHx  %finds rsch at each point
        rstar = x(k);
        if rstar < -1.0e3   %less than -1000, very close to horizon (when m=1)
              error = 'Tortoise coordinate too small'; %#ok<NASGU> 
              pause
        end
        if rstar > twom       %not near horizon
           rsch = rstar;        
           rschm2m = rsch - twom;

         else                %near horizon
          rschm2m = twom*exp((rstar-twom)/2*twom);
          rsch = rschm2m + twom;

        end
        for i=1:100    %Newtons method
          rsn = rschm2m + twom + 2*twom*log(rschm2m/twom) - ...
              (twom^2/rschm2m);
          if abs(rsn-rstar) < 1.0e-25  %approximation is good enough
             break
          end
          drsdr = (rschm2m + twom)^2/(rschm2m^2);
          rtempm2m = rschm2m + (rstar-rsn)/drsdr;
          if rtempm2m < 0
            rtempm2m = 0.1*rschm2m;
          end
          rschm2m  = rtempm2m;
          rsch     = rschm2m + twom;
        end
       Rrsn(k) = rsch;
       Rschm2m(k) = rschm2m;
       
      end