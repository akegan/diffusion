function [ out ] = u_fixedBC(x,t,L,D,A,B)
%1D Concentration, as a function of x,t, L (length of gel), D (diffusion
%coefficient, A (left side boundary concentration), B (right side boundary
%concentration)
%   Detailed explanation goes here
syms n
lambda_n=(n*pi/L)^2;
sum_n=A/(2*pi*n)*exp(-lambda_n*D*t)*sin(n*pi/L*x);
u_fixedBC=symsum(sum_n,n,1,Inf)+x/L*(B-A)+A;
out=u_fixedBC;


end

