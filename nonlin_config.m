function [ options ] = nonlin_config( varargin )
%Set the config of ode solver
%Input:
%       
%Output:
%   options.
%       MaxIteration: Max number of iterations
%       Tol: tolerence
%       Method: method used in solving nonlinear systems
%           NE:mulit-variable newton method
%           Broyden: Broyden method
%       LinMethod: method used in solving linear systems when solving non

%Error trapping:
if(mod(nargin,2)~=0) error('Wrong input variable numbers!');end
%Set default config:
options.MaxIteration=100;
options.Tol=1e-6;
options.Method='Broyden';
options.LinMethod='Gauss';

%Change default config:
for i=1:2:nargin
   switch varargin{i}
       case 'MaxIteration'
           options.MaxIteration=varargin{i+1};
       case 'Tol'
           options.Tol=varargin{i+1};
       case 'Method'
           switch varargin{i+1}
               case 'NE'; options.Method='NE';
               case 'Broyden'; options.Method='Broyden';
               otherwise; error('No such nonlin_sovle method!');
           end
       case'LinMethod'
           switch varargin{i+1}
               case 'Jacobi'; options.LinMethod='Jacobi';
               case 'GS'; options.LinMethod='GS';
               case 'CG'; options.LinMethod='CG';
               case 'Gauss';options.LinMethod='Gauss';
               otherwise; error('No such lin_solve method!');
           end
       otherwise
           error('Input config not found!');
   end
end

end

