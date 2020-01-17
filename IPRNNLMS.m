function [s11]=IPRNNLMS(xn1,dn,w0,w1)  %step-size scalers1
L=30;%taps
N=50000;% iterations
x=xn1;
w=ones(L,N);
w(:,L-1)=w1;
e=zeros(1,N);
beta=3;

u=0.05;
c=ones(L,1);
for n=L:N
    xn=x(n:-1:(n-L+1));
    xn=xn.';
    d=xn.'*xn;
    y=(w(:,n-1))'*xn;
    e(n)=dn(n)-y;
    for o=1:L
        c(o)=w(o,n-1)/(abs(w(o,n-1))+0.01);
    end
    Dc=diag(c);
    s=4*exp(-1*beta*(e(n)^2/d))/((1+exp(-1*beta*(e(n)^2/d)))^2);
   
    w(:,n)=w(:,n-1)+u*s*Dc*e(n)*xn/d;
    
end
s11=zeros(1,N);
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
    s11(j)=10*log10(sum/sum1);
end
s11=s11(L:N);
plot(s11)
hold on