%function [Freq, DR] = ERA(N,n,dt,f)

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

r = round(N/2) - 1;
c1 = f(1:r);
r1 = f(r:N-2);
H0 = hankel(c1,r1);
c2 = f(2:r+1);
r2 = f(r+1:N-1);
H1 = hankel(c2,r2);
[U,S,V] = svd(H0);
Sn = S(1:n,1:n);
Un = U(:,1:n);
Vnt = V(:,1:n);
A = (Sn^(-1/2))*(Un')*(H1)*(Vnt)*(Sn^(-1/2));
E = eig(A);
lambda = log(E)/dt;
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
