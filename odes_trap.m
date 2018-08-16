function [ t,y ] = odes_trap( varargin )
%Solve ode set IVP by Trapezoid method by WangXiao 8/16/2018
%The IVP should be in the following form:
%   yk'=fk(t,y1,y2,...,yn)
%   y(a)=(y1,y2,...,yn)
%   t=[a,b]
%Input:
%(f,inter,y0) or (f,inter,y0,n)
%       f:the derivative of Vector function (in the form cell)
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
n=length(f);
if(n~=length(y0))
    error('Error! demision of derivative and initial value not equal!');
end

%Trapzoid method starts:
N=floor((inter(2)-inter(1))/h+1);
t=linspace(inter(1),inter(2),N);
y=zeros(n,N);
y(:,1)=y0;
k=2;
while(k<=N)
   y(:,k)=transpose(odes_trapstep(f,t(k-1),transpose(y(:,k-1)),h));
   k=k+1;
end

end