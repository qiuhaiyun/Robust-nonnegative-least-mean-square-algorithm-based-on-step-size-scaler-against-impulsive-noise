function [s8]=RNNLMS(xn1,dn,w0,w1)  %step-size scalers1
L=30;%taps
N=50000;% iteration
x=xn1;
w=ones(L,N);
w(:,L-1)=w1;
e=zeros(1,N);
beta=3;

u=0.2;
for n=L:N
    xn=x(n:-1:(n-L+1));
    xn=xn.';
    Dx=diag(xn);
    c=xn.'*xn;
    y=(w(:,n-1))'*xn;
    e(n)=dn(n)-y;
    s=4*exp(-1*beta*(e(n)^2/c))/((1+exp(-1*beta*(e(n)^2/c)))^2);
    
    w(:,n)=w(:,n-1)+u*s*Dx*e(n)*w(:,n-1)/c;
    
end
s8=zeros(1,N);
sum1=0;
for k=1:L
    sum1=sum1+w0(k)^2;
end
for j=L:N
    q=w0-w(:,j);
    sum=0;
    for i=1:L
        sum=sum+q(i)^2;
    end
    s8(j)=10*log10(sum/sum1);
end
s8=s8(L:N);
plot(s8)
hold on