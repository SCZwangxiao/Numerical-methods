function [ options ] = odeconfig( varargin )
%Set the config of ode solver
%Input:
%       
%Output:
%   options.
%       MaxStep: Max step-size for vari-step method/step-size for fixed-step method
%       RelTol: relative tolerence
%       FixStepMethod: method used in fixed step method
%           Euler: Euler  method (1st order)
%           Trap: Trapezoid method £¨2nd order£©
%           RK4: Runge-Kutta method (4th order)

%Error trapping:
if(mod(nargin,2)~=0) error('Wrong input variable numbers!');end
%Set default config:
options.MaxStep=0.1;
options.RelTol=1e-6;
options.FixStepMethod='RK4';

%Change default config:
for i=1:2:nargin
   switch varargin{i}
       case 'MaxStep'
           options.MaxStep=varargin{i+1};
       case 'RelTol'
           options.RelTol=varargin{i+1};
       case 'FixStepMethod'
           switch varargin{i+1}
               case 'Euler'; options.FixStepMethod='Euler';
               case 'Trap'; options.FixStepMethod='Trap';
               case 'RK4'; options.FixStepMethod='RK4';
               otherwise; error('No such fixed-step method!');
           end
       otherwise
           error('Input config not found!');
   end
end

end

