function [dn,w0,w1,xn1]=data4()
L=30;%taps
w1=rand(L,1);   
w2=ones(L,1);
N=50000;%iterations
for i=1:20
    w2(i)=1-0.05*i;
end
for j=21:25
    w2(j)=0;
end
for p=26:30
    w2(p)=-0.01*(p-25);
end
w0=zeros(L,1);
w0(1:20)=w2(1:20);
xn1=ones(1,N);%input
tao=0;%¦Ó=0.5or0
w=normrnd(0,(1-tao^2)^0.5,1,N);
for i=2:N
    xn1(i)=tao*xn1(i-1)+w(i-1);
end
k=zeros(size(xn1));
k(1:L)=xn1(1:L);
for m=L:N
    x=xn1(m:-1:(m-L+1));
    x=x.';
    k(m)=w2.'*x;
end
m=binornd(1,0.1,1,N);%p=0.5
dn1=awgn(k,10,'measured');%SNR=10

%impulsive noise
sum=0;
for i=1:N
    c=mean(dn1);
    sum=sum+(dn1(i)-c)^2;
end
sum=sum/N;
b=normrnd(0,1000*sum,1,N);

dn=dn1+b.*m;