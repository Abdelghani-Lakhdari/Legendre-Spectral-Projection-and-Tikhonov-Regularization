clc;clear all;close all;format long
N=input('N=');N=N+1;epsilon=input('epsilon=');alpha=input('alpha=');
Ga=1/2;g=@(s)2*(1-s)+4/3*(s^(3/2))*(1-s)^0.5-4/3*(1-s)^2;
f=@(t)t;K=@(s,t)(1-s)^0.5;a=0;b=1;V=((b-a)/2)^(1-Ga);M=N;
Ps=@(x) (b-a)/2*x+(b+a)/2;Ts=@(t,x) (1+t)/2*x+(t-1)/2;Tr=@(t,x) (1-t)/2*x+(t+1)/2;
[E,WW] =  GaussJacobi(M,0,0);[T1,W1] = GaussJacobi(M,0,-Ga);[T2,W2] = GaussJacobi(M,-Ga,0);
for j=1:M
for i=1:N
S0=0;S1=0;
for k=1:M
A1(i,j)=S0+V*((1+E(j))/2)^(1-Ga)*W1(k)*K(Ps(E(j)),Ps(Ts(E(j),T1(k))))*sqrt((2*i-1)/2)*Legn(i-1,Ts(E(j),T1(k)));
S0=A1(i,j);
A2(i,j)=S1+V*((1-E(j))/2)^(1-Ga)*W2(k)*K(Ps(E(j)),Ps(Tr(E(j),T2(k))))*sqrt((2*i-1)/2)*Legn(i-1,Tr(E(j),T2(k)));
S1=A2(i,j);
end       
end
end
PSI=A1+A2;
for j=1:N
for i=1:N
s0=0;
for k=1:M
A(i,j)=s0+WW(k)*PSI(i,k)*PSI(j,k);s0=A(i,j);
end       
end
end

for k=1:M
G(k)=g(Ps(E(k)));
end
RAN=randn(size(G));GG1=G+norm(G)*epsilon*RAN/norm(RAN);
for j=1:N
s=0;
for k=1:M
GG(j)=s+WW(k)*GG1(k)*PSI(j,k);s=GG(j);
end
end
X=-1:0.02:1;[~,MM]=size(X);
for j=1:MM
for i=1:N
PSI2(i,j)=sqrt((2*i-1)/2)*Legn(i-1,X(j));
end       
end 
X=Ps(X);delta=norm(G-GG1);coef=(alpha*eye(size(A))+A)\GG';
for j=1:MM
s=0;
for k=1:N
Fapp(j)=s+coef(k)*PSI2(k,j);s=Fapp(j);
end
end
for j=1:MM
Fex(j)=f(X(j));
end
err=abs(Fex-Fapp);
Er=max(err)
Rel=norm(err)/norm(Fex)
figure(1)
hold on,
subplot(2,1,1)
plot(X,Fex,'--rs',X,Fapp,'b-', 'LineWidth',1,'MarkerSize',2),
    grid on,
    legend('Excact solution','Approximate solution'),
    subplot(2,1,2)
    plot(X,err,'r','LineWidth',1),
    legend(strcat('\delta=',num2str(epsilon),', n =',num2str(N-1),', \alpha =',num2str(alpha))), 
    title('Error = |Excact - Approximate |'),
    grid on,