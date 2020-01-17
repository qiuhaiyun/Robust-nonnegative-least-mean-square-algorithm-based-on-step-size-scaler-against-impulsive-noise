function [dn,w0,w1,xn1]=data2()
L=10;%taps
w1=rand(L,1);   
w0=[0.8 0.6 0.5 0.4 0.3 0.2 0.1 0 0 0].';
w2=[0.8 0.6 0.5 0.4 0.3 0.2 0.1 -0.1 -0.3 -0.6].';
N=50000;% iteration
xn1=normrnd(0,1,1,N);%input
k=zeros(size(xn1));
k(1:L)=xn1(1:L);
for m=L:N
    x=xn1(m:-1:(m-L+1));
    x=x.';
    k(m)=w2.'*x;
end
m=binornd(1,0.5,1,N);%p=0.5
%noise=2*rand(1,N)-1;
dn1=awgn(k,10,'measured');%dn1=k+noise;%SNR=10

%impulsive noise
sum=0;
for i=1:N
    c=mean(dn1);
    sum=sum+(dn1(i)-c)^2;
end
sum=sum/N;
b=normrnd(0,1000*sum,1,N);
dn=dn1+b.*m;