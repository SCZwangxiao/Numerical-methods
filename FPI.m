function [varargout] = FPI(varargin)
%{
Fixed point iteration method by Wang Xiao 8.12.2018
Input:
(g,p0,Tol) or (g,p0,Tol,N)
        g:The function of g(x)=x
        p0:guess of fixed point
        Tol:tolerence of absolute error (optional)
        N:the maximun number of iteration (optional)
Output:
(pc) or (pc,iter)
        pc:The approximation of zero
        iter:Struct shows the process of the iteration
            iter.zero approximation zero in each iteration
            iter.N times of iteration
            iter.FP approximation fixed point in each iteration
            iter.Xerror: xerror |pn-pn-1| in each iteration
%}

%defaultsettings
Tol=1e-8;
N=100;

%process input
if(~(nargout==1||nargout==2))
    error('Error! output variable numbers incorrect!');
end
if(nargin==2)
    g=varargin{1};
    p0=varargin{2};
elseif(nargin==3)
    g=varargin{1};
    p0=varargin{2};
    Tol=varargin{3};
elseif(nargin==4)
    g=varargin{1};
    p0=varargin{2};
    Tol=varargin{3};
    N=varargin{4};
else
    error('Error! input variable numbers incorrect!');
end

%Fixed point method starts:
xerror=zeros(1,N);
p=zeros(1,N);

p(1)=g(p0);
xerror(1)=abs(p(1)-p0);
i=2;
while(i<=N)
    p(i)=g(p(i-1));
    p(i)
    xerror(i)=abs(p(i)-p(i-1));
    if(xerror(i)<Tol)
       break;
    end
    i=i+1;
end

if(i==N+1)
    error('Error! Maximum numbers of iteration exceed!');
else
   pc=p(i); 
   Nmax=i;
end

if(nargout==1)
    varargout{1}=pc;
elseif(nargout==2)
    varargout{1}=pc;
    iter.zero=p(1:Nmax);
    iter.N=Nmax;
    iter.FP=p(1:Nmax);
    iter.Xerror=xerror(1:Nmax);
    varargout{2}=iter;
else
    error('Error! output variable numbers incorrect!');
end

end

