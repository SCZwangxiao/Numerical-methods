function [ options ] = lineqconfig( varargin )
%Set the config of linear equation solver
%Input:
%       
%Output:
%   options.
%       MaxIteration: Max Iterations
%       Tol: tolerence
%       Method: method used
%           Jacobi: Jacobi method
%           GS: Gauss-Seidel method
%           CG: Conjugate Gradient
%       SOR: set the SOR parameter w

%Error trapping:
if(mod(nargin,2)~=0);error('Wrong input variable numbers!');end
%Set default config:
options.MaxIteration=100;
options.Tol=1e-6;
options.Method='GS';
options.SOR=1;

%Change default config:
for i=1:2:nargin
   switch varargin{i}
       case 'MaxIteration'
           options.MaxIteration=varargin{i+1};
       case 'Tol'
           options.Tol=varargin{i+1};
       case 'Method'
           switch varargin{i+1}
               case 'Jacobi'; options.Method='Jacobi';
               case 'GS'; options.Method='GS';
               case 'CG'; options.Method='CG';
               otherwise; error('No such iteration method!');
           end
       case 'SOR'
           options.SOR=varargin{i+1};
           if(options.SOR<0||options.SOR>2); error('w in SOR method must be between[0 2]!');end
       otherwise
           error('Input config not found!');
   end
end

end



