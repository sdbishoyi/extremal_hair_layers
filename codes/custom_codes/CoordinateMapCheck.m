clear
x=[-100:1:100]; M=1; a=0.01;
 
 rm = M - sqrt(M^2 - a^2)  
 rp = M + sqrt(M^2 - a^2) 
 
 rkerr(:,1)=x';
 rkerr(:,2)=RstarToR(x,M)';  %%% scott's code
 rkerr(:,3)=RstarToRwithSpin(x,M,a)';   %%% my code using NR
 rkerr(:,4)=RstarToRwithSpin2(x,M,a)';  %%% my code using fzero
 
%  RstarToRwithSpin2(x,M,a)';