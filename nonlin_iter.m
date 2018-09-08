function [ varargout ] = nonlin_iter( varargin )
%{
Solving nonlinear systems by itertion method
By Wang Xiao 9/7/2018
The nonlinear problem should be in the form:
    F(x1,x2,...,xn)=0
Input:
        F: (Vector valued) function to solve [function handle] [column vector]
        DF/B: the Jacobian or approximate one of F
        x0: initial value of (Vector valued) function [column vector]
        [options:a structure array thats show the config of the ode solver
Output:
        xp: estimated zero [column vector]
        iter: a structure that contains results
            x: 
            iterations:
            
%}

%Input variables:
if(nargin==3)
    F=varargin{1};
    B=varargin{2};
    x0=varargin{3};
    options=odeconfig();
elseif(nargin==4)
    F=varargin{1};
    x0=varargin{3};
    options=varargin{4};
    if(strcmp(options.Method,'NE')); DF=varargin{2};
    else B=varargin{2}; end
else
    error('Wrong input variable numbers!')
end

%Error trapping:


%Fixed-step method starts:
%Set parameters
n=length(x0);   % number of equations
N=options.MaxIteration;
Tol=options.Tol;
k=1;    % iterations that finished -1

%Initialize
x=zeros(n,N+1);
x(:,1)=x0;

%Iteration starts
switch options.Method
    case 'NE'
        switch options.LinMethod
            case 'Gauss' 
                while(k<=N)
                    s=-DF(x(:,k))\F(x(:,k));
                    x(:,k+1)=x(:,k)+s;
                    k=k+1;
                    if(norm(s, Inf)<Tol); break; end
                end    
            otherwise
                linoptions=lineqconfig('Method',options.LinMethod,'Tol',Tol);
                while(k<=N)
                    s=lineq_iter( DF(x(:,k)), -F(x(:,k)), ones(n),linoptions);
                    x(:,k+1)=x(:,k)+s;
                    k=k+1;
                    if(norm(s, Inf)<Tol); break; end
                end
        end
        xp=x(:,k);
    case 'Broyden'     
        while(k<=N)
            s=-B*F(x(:,k));
            x(:,k+1)=x(:,k)+s;
            if(norm(s,Inf)<Tol); k=k+1;break; end
            delta=x(:,k+1)-x(:,k);
            DELTA=F(x(:,k+1))-F(x(:,k));
            B=B+(delta-B*DELTA)*transpose(delta)*B / (transpose(delta)*B*DELTA);        
            k=k+1;
        end
        xp=x(:,k);
end

if(k==N+1)
    disp('The method ends at:');
    disp(xp);
    error('Max Iterations exceeds, solution not found');
end

if(nargout==1)
    varargout{1}=xp;
elseif(nargout==2)
    varargout{1}=xp;
    iter.x=x(:,2:k+1);
    iter.iterations=k-1;
    varargout{2}=iter;
else
end
    
end

