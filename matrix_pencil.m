%function [Freq, DR] = matrix_pencil(N,n,dt,y)

clc
clear all
load kundur_all_methods
aux2 = find(t==6);
y = Pe3_1(aux2:end);
t2 = t(aux2:end);
plot(t2,y);
N = 200%length(f);
n = 4;
dt = 1/60;

%Construimos la matriz hankel H0
L1 = ceil(1/3 * N);
L2 = floor(2/3 * N);
L = ceil((L1 + L2) / 2); % obtenemos el parametro de pencil
for i=1:(N-L)
    Y(i,:)=y(i:(i+L));
end
[U,S,V] = svd(Y); %Descomponemos la matriz hankel

Vnew=V(:,1:n);
a=size(Vnew,1);
Vs1=Vnew(1:(a-1),:);
Vs2=Vnew(2:end,:);
Y1=Vs1'*Vs1;   
Y2=Vs2'*Vs1;
Y1=Y1';
E = eig((inv(Y1))*Y2);
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