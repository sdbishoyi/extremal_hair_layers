%%This is the final code to be used for converting rstar to r 
%%written by Som on 12-22-2023. Used the suggestion from Gaurav to write the algorithm from right to left
%%And Scott suggested to use isnan function to check for NaNs and assign rp at that point


%%%% Katie pointed out on 20/01 meeting that the file name did not have _Newton but the code works perfectly in notebooks
%%%% one more thing that there was an end missing for the first for loop
%%%% I do not have any idea but looking for it--- as of 24/01/22

function [rkerr] = RstarToRwithCharge_Newton(rs,m)

%     if length(rs)>=1
    
                rkerr = zeros(size(rs));
                rp = m;
                len = size(rs,1)*size(rs,2);
                r=rs(len);
            
                for k=1:len
                    rstar=rs(len-k+1);
                    
               
                    for i=1:100   %Newtons method
                              funct     = r + 2*m*log(abs(r-rp)/rp) ...
                                            - (m^2./(r - rp)) ...
                                            - rstar;
                              fprime    = (1 - m./r).^(-2)  ;
                                            
                              r = r - funct/fprime;
                              
                        if isnan(r)==1
                                    r=rp+1e-15;
                        end
                    end
                       
                   
                    %if rstar<-50
                    %   r=rp;
                    %end
                    if r<rp
                        r=rp+1e-14;
                    end    
                    
                   
                    rkerr(len-k+1)=r;
                end

                rkerr(1)=rp;

%     else
%                 rm = m - sqrt(m^2 - a^2) ; 
%                 rp = m + sqrt(m^2 - a^2) ;
%                 r=rs;
%                 rkerr=r - (2*m*rp/(rp-rm))*log(abs(r-rp)/2*m) ...
%                                             + (2*m*rm/(rp-rm))*log(abs(r-rm)/2*m);
% 
%     end

end
    
  
