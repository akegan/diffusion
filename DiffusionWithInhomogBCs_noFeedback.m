clear all;
A=5;
B=0;
D=.1;
L=1;
tfinal=4;
Nmax=500;
dx=.01;
dt=.01;
x=0:dx:L;
t=0:dt:tfinal;
n=1:1:Nmax;
[X,T,N]=ndgrid(x,t,n);
[x2,t2]=ndgrid(x,t);
sum_n=-sum(2*(A+B*cos(pi*N))./(pi*N).*exp(-(N*pi/L).^2*D.*T).*sin(N*pi/L.*X),3)+x2/L.*(B-A)+A;
%test=exp(-(n*pi/L).^2*D.*t);%.*
plot(x,sum_n(:,1:50),'b')
sum_n(100,2)
