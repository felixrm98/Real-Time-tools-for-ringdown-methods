%function [Freq, DR] = prony(N,n,dt,f)
clc
clear all
load kundur_all_methods
aux2 = find(t==6);
f = Pe3_1(aux2:end);
t2 = t(aux2:end);
plot(t2,f);
N = 200%length(f);
n = 4;
dt = 1/60;


%f = f - mean(f);
b = f(n+1:N);
b1 = f(n:N-1);
R1 = f(n:-1:1);
T = toeplitz(b1,R1);
a = pinv(T)*b;
p = [1;(a(1:n)*-1)];
z = roots(p);
lambda = log(z)/dt;
sigma = -real(lambda);
omega = imag(lambda);
fre1 = omega/(2*pi);
damp = sigma;
aux = [fre1 damp];
aux = sortrows(aux);
fre1 = aux(:,1);
damp = aux(:,2);
Freq = fre1((n/2)+1:n)
dr1 = damp((n/2)+1:n);
DR =(dr1./sqrt((dr1.^2)+((Freq*2*pi).^2)))*100
%end
