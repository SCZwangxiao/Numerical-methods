function [ t,y ] = ode_fixstep( varargin )
%Solving IVP of odes using fixed-step method.
%By Wang Xiao 8/29/2018
%The IVP should be in the following form:
%	yk'=fk(t,y1,y2,...,yn) ,k=1,2...n
%   y(a)=(y1,y2,...,yn) [column vector]
%   t=[a,b]
%Input:
%       odefun: (Vector valued) function to solve [function handle] [column vector]
%           the function must be defined as f=odefun(t,y)
%       tspan: interval of the integration, tspan=[a,b]
%       y0: initial value of (Vector valued) function [column vector]
%       [options:a structure array thats show the config of the ode solver
%Output:
%       t: evalution points [column vector]
%       y: solution of function [column vector]

%Input variables:
if(nargin==3)
    odefun=varargin{1};
    tspan=varargin{2};
    y0=varargin{3};
    options=odeconfig();
elseif(nargin==4)
    odefun=varargin{1};
    tspan=varargin{2};
    y0=varargin{3};
    options=varargin{4};
else
    error('Wrong input variable numbers!')
end

%Error trapping:
if(length(tspan)~=2); error('wrong input of tspan!');end

%Fixed-step method starts:
%Set parameters
n=length(y0);   % number of odes
h=options.MaxStep;
N=floor((tspan(2)-tspan(1))/h)+1; % number of iterations
h=(tspan(2)-tspan(1))/(N-1);
if(h~=options.MaxStep); fprintf('Step-size has bee changed to %.8f\n',h);end
%Initialize t,y
t=linspace(tspan(1),tspan(2),N);
y=zeros(n,N);
y(:,1)=y0;
%Iteration starts
switch options.FixStepMethod
    case 'Euler'
        for k=1:N-1
        y(:,k+1)=y(:,k)+h*odefun(t(k),y(:,k));
        end
    case 'Trap'
        for k=1:N-1
        s1=odefun(t(k),y(:,k));
        s2=odefun(t(k+1),y(:,k)+h*s1);
        y(:,k+1)=y(:,k)+h*(s1+s2)/2;
        end
    case 'RK4'
        for k=1:N-1
        s1=odefun(t(k),y(:,k));
        s2=odefun(t(k)+h/2,y(:,k)+h*s1/2);
        s3=odefun(t(k)+h/2,y(:,k)+h*s2/2);
        s4=odefun(t(k+1),y(:,k)+h*s3);
        y(:,k+1)=y(:,k)+h*(s1+2*s2+2*s3+s4)/6;  
        end
end

end

