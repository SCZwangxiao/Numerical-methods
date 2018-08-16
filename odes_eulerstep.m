function [ yb ] = odes_eulerstep( f,ta,ya,h )
%One step in the Euler method by WangXiao 8/16/2018
%The method should be in the following form:
%   yb=ya+h*{f(ta,ya)+f[ta+h,ya+h*f(ta,ya)]}/2
%Input:
%(f,ta,ya,h)
%       f:the derivative of Vector function
%       ta:Value of independent variable in the current step
%       ya:Vector of value function in the current step (n demision)
%       h:the length of steps
%Output:
%       yb:Certain component of Vector function in the next step
n=length(ya);
if(length(f)~=n)s
    error('Error! demision of derivative and function not equal!');
end
yb=zeros(1,n);
for k=1:n
ff(k)=f{k}(ta,ya);
end
i=1;
while(i<=n)
    yb(i)=ya(i)+h*ff(i);
    i=i+1;
end

end

