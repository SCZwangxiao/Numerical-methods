function [ varargout ] = lineq_iter( varargin )
%{
Solving linear systems using iteration method
By Wang Xiao  8/30/2018 
The problem should be in the form Ax=b 
    A: n*n matrix
    b: n dimension column vector
Input:
        A: the coefficient matrix
        b: the vector [column vector]
        x0: the initial guess [column vector]
        [options:a structure array thats show the config of the lineq solver
Output:
        xp: the approximate solution [column vector]
        iter: the struct that contains:
            iterations: number of iterations
            x: a sequence of the previous results
%}

%Input variables:
if(nargin==3)
    A=varargin{1};
    b=varargin{2};
    x0=varargin{3};
    options=lineqconfig();
elseif(nargin==4)
    A=varargin{1};
    b=varargin{2};
    x0=varargin{3};
    options=varargin{4};
else
    error('Wrong input variable numbers!');
end
    
%Error trapping:
[numA_row,numA_col]=size(A);
if(numA_row~=numA_col);error('Input coefficient matrix is not square!');end
[numb_row,numb_col]=size(b);
if(numb_row~=numA_row||numb_col~=1);error('Input b is not in right dimension');end
[numx0_row,numx0_col]=size(x0);
if(numx0_row~=numA_row||numx0_col~=1);error('Input x0 is not in right dimension');end

%Set parameters
N=options.MaxIteration;
n=numA_row;
Tol=options.Tol;
SORw=options.SOR;
%Initialize 
xp=x0;
k=0; % the number of iteration that already finished
x=zeros(n,N);

%Iteration starts:
switch options.Method
%Jacobi method
    case 'Jacobi'
        %Judge the convergence
        %{
        for i=1:n
           sum=0;
           for j=[1:i-1,i+1:n];sum=sum+abs(A(i,j));end
           if( abs(A(i,i))<=sum ); disp('The coefficient matrix is not strictly-diagonally-dominant,the method may not converge!');break;end
        end
        %}

        %Iteration
        k=0;
        while(k<N)
            D=diag(A);
            r=A-diag(D); %(matirx L+U)
            xp=(b-r*x0)./D;
            if( norm((xp-x0),Inf)<Tol);break;
            else x0=xp;end
            k=k+1;
            x(:,k)=xp;
        end
        if(k==N)
            disp('The method ends at:');
            disp(xp);
            error('Max Iterations exceeds, solution not found');
        end;
%Gauss-Seidel method
    case 'GS'
        %Judge the convergence
        for i=1:n
           sum=0;
           for j=[1:i-1,i+1:n];sum=sum+abs(A(i,j));end
           if( abs(A(i,i))<=sum ); disp('The coefficient matrix is not strictly-diagonally-dominant,the method may not converge!');break;end
        end
        %Iteration
        k=0;
        if(SORw==1)
            while(k<N)
                for i=1:n
                    sum=0;
                    for j=1:i-1; sum=sum-A(i,j)*xp(j);end
                    for j=i+1:n; sum=sum-A(i,j)*x0(j);end
                    xp(i)=(sum+b(i))/A(i,i);
                end
                if( norm((xp-x0),Inf)<Tol ); break;
                else x0=xp;end
                k=k+1; 
                x(:,k)=xp;
            end
        else
            SORw_=(1-SORw);
            while(k<N)
                for i=1:n
                    sum=0;
                    for j=1:i-1; sum=sum-A(i,j)*xp(j);end
                    for j=i+1:n; sum=sum-A(i,j)*x0(j);end
                    xp(i)=SORw_*x0(i)+(sum+b(i))*SORw/A(i,i);
                end
                if( norm((xp-x0),Inf)<Tol ); break;
                else x0=xp;end
                k=k+1; 
                x(:,k)=xp;
            end
        end
        
        if(k==N)
            disp('The method ends at:');
            disp(xp);
            error('Max Iterations exceeds, solution not found');
        end
%Conjugate Gradient method
    case 'CG'
        d=zeros(n,N);
        r=zeros(n,N);
        alpha=zeros(1,N);
        beta=zeros(1,N);
        %initialize: 
        r0=b-A*x0;
        d0=r0;
        
        %1st iteration
        alpha0=transpose(r0)*r0/( transpose(d0)*A*d0 );
        x(:,k+1)=x0+alpha0*d0;
        r(:,k+1)=r0-alpha0*A*d0;
        beta0=transpose(r(:,k+1))*r(:,k+1)/ ( transpose(r0)*r0 );
        d(:,k+1)=r(:,k+1)+beta0*d0;
        k=k+1;
        %2-n iteration
        while(k<N)
           alpha(k)=transpose(r(:,k))*r(:,k)/( transpose(d(:,k))*A*d(:,k) );
           x(:,k+1)=x(:,k)+alpha(k)*d(:,k);
           r(:,k+1)=r(:,k)-alpha(k)*A*d(:,k);
           beta(k)=transpose(r(:,k+1))*r(:,k+1)/ ( transpose(r(:,k))*r(:,k) );
           d(:,k+1)=r(:,k+1)+beta(k)*d(:,k);
           k=k+1;
           if( norm((r(:,k)-r(:,k-1)),Inf)<Tol ); break; end
        end
        xp=x(:,k);
        if(k==N); 
            disp('The method ends at:');
            disp(xp);
            error('Max Iterations exceeds, solution not found');
        end
end

if(nargout==1)
    varargout{1}=xp;
elseif(nargout==2)
    varargout{1}=xp;
    iter.x=x(:,1:k);
    iter.iterations=k;
    varargout{2}=iter;
end

