function [varargout] = bisect(varargin)
%{
Bisection method by Wang Xiao 8.10.2018
Input:
(f,a,b) or (f,a,b,Tol) or (f,a,b,Tol,N)
        f:The function of f(x)=0
        a:left border of guess
        b:right border of guess
        Tol:tolerence of absolute error (optional)
        N:the maximun number of iteration (optional)
Output:
(xc) or (xc,iter)
        xc:The approximation of zero
        iter:Struct shows the process of the iteration
            iter.N times of iteration
            iter.zero approximation zero in each iteration
            iter.Xerror: xerror |xn-xn-1| in each iteration
            iter.Xerror: yerror |f(xn)| in each iteration
%}

%defaultsettings
Tol=1e-8;
N=100;

%process input
if(~(nargout==1||nargout==2))
    error('Error! output variable numbers incorrect!');
end
if(nargin==3)
    f=varargin{1};
    a=varargin{2};
    b=varargin{3};
elseif(nargin==4)
    f=varargin{1};
    a=varargin{2};
    b=varargin{3};
    Tol=varargin{4};
elseif(nargin==5)
    f=varargin{1};
    a=varargin{2};
    b=varargin{3};
    Tol=varargin{4};
    N=varargin{5};
else
    error('Error! input variable numbers incorrect!');
end

%Biection method starts:
if(f(a)*f(b)>0)
   error('Error! f(a)*f(b)>0');
end

FP=zeros(1,N);
xerror=zeros(1,N);
c=zeros(1,N);

i=1;
while(i<=N)
   c(i)=(a+b)/2;
   FP(i)=f(c(i));
   xerror(i)=abs(b-a)/2;
   if(FP(i)==0 || xerror(i)<Tol)
       break;
   end
   
   if(FP(i)*f(a)<0)
       b=c(i);
   elseif(FP(i)*f(b)<0)
       a=c(i);
   end
   i=i+1;
end
if(i==N+1)
    error('Error! Maximum numbers of iteration exceed!');
else
   xc=c(i); 
   Nmax=i;
end

if(nargout==1)
    varargout{1}=xc;
elseif(nargout==2)
    varargout{1}=xc;
    iter.N=Nmax;
    iter.zero=c(1:Nmax);
    iter.Xerror=xerror(1:Nmax);
    iter.Yerror=FP(1:Nmax);
    varargout{2}=iter;
else
    error('Error! output variable numbers incorrect!');
end

end