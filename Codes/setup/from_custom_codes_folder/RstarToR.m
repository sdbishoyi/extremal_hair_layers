function [Rsch,Rschm2m] = RstarToR(x,m)

%Uses newton's method to turn tortoise coordinates to standard schwarschild
%coordinates
%rstar = read in rstar coordinates, rsch = schwarzchild r coordinate..but
%computed iteratively using newtons method
%rschm2m + 2M ~ rsch... so rschm2m near 0 at horizon instead of 2m,
%we actually update rshm2m and add 2m to this number
%rsn = rstar (newton)

      Rsch = zeros(size(x));
      Rschm2m = zeros(size(x));
      LENGTHx = size(x,1)*size(x,2);
      twom   = 2*m;
      for k = 1:LENGTHx  %finds rsch at each point
        rstar = x(k);
        if rstar < -1.0e3   %less than -1000, very close to horizon (when m=1)
              error = 'Tortoise coordinate too small'
              pause
        end
        if rstar > 4.0*m   %not near horizon
          rsch    = rstar;
          rschm2m = rsch - twom;
        else                %near horizon
          rschm2m = twom*exp(rstar/twom - 1);
          rsch    = rschm2m + twom;
        end
        for i=1:100    %Newtons method
          rsn   = rschm2m + twom*(1+log(rschm2m/twom));
          if abs(rsn-rstar) < 1.0e-25  %approximation is good enough
             break
          end
          drsdr = (rschm2m+twom)/rschm2m;
          rtempm2m = rschm2m + (rstar-rsn)/drsdr;
          if rtempm2m < 0
            rtempm2m = 0.5*rschm2m;
          end
          rschm2m  = rtempm2m;
          rsch     = rschm2m + twom;
        end
       Rsch(k) = rsch;
       Rschm2m(k) = rschm2m;
       
      end
