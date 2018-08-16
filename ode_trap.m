function [ t,y ] = ode_trap( varargin )
%Solve ode IVP by Trapezoid method by WangXiao 8/16/2018
%The IVP should be in the following form:
%   y'=f(t,y)
%   y(a)=y0
%   t=[a,b]
%Input:
%(f,inter,y0) or (f,inter,y0,n)
%       f:the derivative of y
%       inter:interval [a,b]
%       y0:initial point
%       h:the length of steps (optional)
%Output:
%       t:
%       y:Solution
%Default settings:
h=0.01;

%Input process
if(nargin==3)
    f=varargin{1};
    inter=varargin{2};
    y0=varargin{3};
    
elseif(nargin==4)
    f=varargin{1};
    inter=varargin{2};
    y0=varargin{3};
    h=varargin{4};
else
    error('Error! input variable numbers incorrect!')
end

%Trapzoid method starts:
n=floor((inter(2)-inter(1))/h+1);
t=linspace(inter(1),inter(2),n);
t=[t,t(n)+h];
y=zeros(1,n);
y(1)=y0;
i=2;
while(i<=n)
   ff=f(t(i-1),y(i-1));
   y(i)=y(i-1)+h*(ff+f(t(i),y(i-1)+h*ff))/2; 
   i=i+1;
end
t=t(1:n);
end