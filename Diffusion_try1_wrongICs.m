clear all;
A=50;
B=0;
D=.01;
L=1;
Nmax=4;
x=0:.001:1;
t=0:1:3;
n=1:1:Nmax;
[X,T,N]=ndgrid(x,t,n);
[x2,t2]=ndgrid(x,t);
sum_n=sum(A./(2*pi*N).*exp(-(N*pi/L).^2*D.*T).*sin(N*pi/L.*X),3)+x2/L.*(B-A)+A;
%test=exp(-(n*pi/L).^2*D.*t);%.*
plot(x,sum_n)