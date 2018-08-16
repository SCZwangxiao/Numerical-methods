function [ t,y ] = ode_euler( varargin )
%Solve ode IVP by Euler method by WangXiao 8/16/2018
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

%Euler method starts:
n=floor((inter(2)-inter(1))/h+1);
t=linspace(inter(1),inter(2),n);
y=zeros(1,n);
y(1)=y0;
i=2;
while(i<=n)
   y(i)=y(i-1)+h*f(t(i-1),y(i-1)); 
   i=i+1;
end
end

