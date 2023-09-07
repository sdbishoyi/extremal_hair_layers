a=0;
b=2;
x=linspace(a,b,5);
% [y,z]=lgwt(3,a,b);
[y,z]=lgwt(length(x),b,a);
% f=;
% z
z=-z;
f=y.^2;
% f=sin(x);

sum(sum(f.*z));

% x, sin(x)
% y, sin(y)
% y
% f=y.^2;
% f=cos(x);
% sum(f.*z);
% sum(sum(f.*z))