n=9;
tspan=[0 1];
bound=[0 exp(1)/3];
h=(tspan(2)-tspan(1)) / (n+1);


e = ones(n,1);
ODE_A=spdiags([e -(h^2+2)*e e],-1:1,n,n);
ODE_b=zeros(n,1);
for i=1:n
   ODE_b(i)=2*h^2*exp(h*i)/3; 
end
ODE_b(1)=ODE_b(1)-bound(1);
ODE_b(n)=ODE_b(n)-bound(2);

w=zeros(n,1);
options=lineqconfig('Method','CG','MaxIteration',2000,'Tol',1e-6);
[xp,iter]=lineq_iter(ODE_A,ODE_b,w,options);
%xp=ODE_A\ODE_b;
t=tspan(1):h:tspan(2);
w=[bound(1);xp;bound(2)];
y=t.*exp(t)/3;
plot(t,w,'o-',t,y,'-')
