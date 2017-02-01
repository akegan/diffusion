clear all;
%%Parameters of the problem
gamma=.25; %partition coefficient 
C_L=1; %Left boundary condition concentration outside gel
D=2; %Diffusion coefficient
V=5; %Volume of right-side well
L=1; %length of gel

%Take into account partition coefficient
A=gamma*C_L; %LHS effective concentration
%D=; %D_effective= 
V=V/gamma; %v_effective=gamma*V


%%Parameters of the calculation
dt=.01; %y step
dx=.01; %x step
Nmax=100; %number of terms to calculate in sum
tfinal=50; %end time

x=0:dx:L;
t=0:dt:tfinal;
n=1:1:Nmax;

u=zeros(size(x,2),size(t,2));
B=0; %initial RHS concentration (also fixed in solution->if I want a different B(0), need to fix)
beta=0; %dB/dt component of an
dx_L=0; %initial condition
B_time=zeros(size(t,2),1);
[X,N]=ndgrid(x,n);
[X2]=ndgrid(x);

%Step over time and calculate u(x,t)
for i=1:size(t,2)
    u_i=A+X2/L*(B-A)+sum((exp(-(N*pi/L).^2*D*t(i)).*(-2*A./(pi*N))+2*(cos(N*pi)./(pi*N)).*beta).*sin(N*pi/L.*X),2); %concentration at a given time step
    dx_L=-(u_i(size(x,2)-1)-u_i(size(x,2)))/dx; %x-derivative of u at x=L boundary
    B=B+dt*(-D/V)*dx_L; %time step RHS boundary
    B_time(i)=B;
    beta=exp(-(N*pi/L).^2*D*t(i)).*beta+dt*(-D/V)*dx_L; %incrementing of dB/dt contribution to An's
    u(:,i)=u_i; %store the time-step data
end
plotstep=1/dt; %plot in 1 second increments
%plot(x,u(:,1:plotstep:size(t,2))) 
    