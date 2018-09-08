function [varargout] = NE(varargin)
%{
Newton method by,Wang Xiao 8.12.2018
Input:
(f,ff,x0) or (f,ff,x0,Tol) or (f,ff,x0,Tol,N)
        f:The function of f(x)=0;
        ff:The derivaive of f(x);
        x0:guess of zero
        Tol:tolerence ofabsolute error(default)
        N:the maximun number of iteration(default)
Output:
        xc:The approximation of zero
        iter:Struct shows the process of the iteration
            iter
            iter.N times of iteration
            iter.Xerror: xerror |xn-xn-1| in each iteration
            iter.Yerror: yerror |f(xn)| approximation  in each iteration
%}

%defaultsettings
Tol=1e-8;
N=1000;

%Process input:
if(~(nargout==1||nargout==2))
    error('Error! output variable numbers incorrect!');
end
if(nargin==3)
    f=varargin{1};
    ff=varargin{2};
    x0=varargin{3};
elseif(nargin==4)
    f=varargin{1};
    ff=varargin{2};
    x0=varargin{3};
    Tol=varargin{4};
elseif(nargin==5)
    f=varargin{1};
    ff=varargin{2};
    x0=varargin{3};
    Tol=varargin{4};
    N=varargin{5};
else
    error('Error! input variable numbers incorrect!');
end

%Fixed point method starts:
xerror=zeros(1,N);
yerror=zeros(1,N);
x=zeros(1,N);

x(1)=x0-f(x0)/ff(x0);
xerror(1)=abs(x(1)-x0);
yerror(1)=f(x(1));
i=2;
while(i<=N)
    x(i)=x(i-1)-yerror(i-1)/ff(x(i-1));
    xerror(i)=abs(x(i)-x(i-1));
    yerror(i)=f(x(i));
    if(xerror(i)<Tol)
       break;
    end
    i=i+1;
end
if(i==N+1)
    error('Error! Maximum numbers of iteration exceed!');
else
   xc=x(i);
   Nmax=i;
end

%Output results:
if(nargout==1)
    varargout{1}=xc;
elseif(nargout==2)
    varargout{1}=xc;
    iter.zero=x(1:Nmax);
    iter.N=Nmax;
    iter.Xerror=xerror(1:Nmax);
    iter.Yerror=yerror(1:Nmax);
    varargout{2}=iter;
else
    error('Error! output variable numbers incorrect!');
end

end

