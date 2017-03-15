clear all

L=1;
%ni=9;
ni=31;
cfstr='ss';
x=0:L/(ni-1):L;n=numel(x);
%%

lsobreh=10;
h=L/lsobreh;b=1;
carga=1;A=b*h;k=5/6;
E=30e6;I=b*h^3/12;G=E/(2*(1+0.3));rho=1;
%%

c=L/sqrt(numel(x));

[xi,xj]=meshgrid(x,x);
AA=g(c,xi,xj);
DAA=dgdx(c,xi,xj);
D2AA=d2gdx2(c,xi,xj);

digits(500)

% 
%equivalente ao dmatriz
% S_total=[(G*A*k)*D2AA,(G*A*k)*DAA ; -G*A*k*DAA,E*I*D2AA-G*A*k*AA];% w, theta
% A_total=[-rho*A*AA,zeros(n,n);zeros(n,n),-rho*I*AA];

S_total(1:n,1:n)=(G*A*k)*D2AA;
S_total(1:n,n+1:2*n)=(G*A*k)*DAA;

S_total(n+1:2*n,1:n)=-G*A*k*DAA;
S_total(n+1:2*n,n+1:2*n)=E*I*D2AA-G*A*k*AA;

A_total(1:n,1:n)=-rho*A*AA;
A_total(1:n,n+1:2*n)=0;

A_total(n+1:2*n,1:n)=0;
A_total(n+1:2*n,n+1:2*n)=-rho*I*AA;


switch cfstr
    
    case {'ss'}
%equacao de fronteira
b=find(x==0 | x==L );
S_total(b,:)=[AA(b,:),zeros(2,n)];
S_total(b+n,:)=[zeros(2,n),E*I*DAA(b,:)];
A_total(b,:)=[zeros(2,n),zeros(2,n)];
A_total(b+n,:)=[zeros(2,n),zeros(2,n)];
    case {'cc'}
b=find(x==0 | x==L );
S_total(b,:)=[AA(b,:),zeros(2,n)];
S_total(b+n,:)=[zeros(2,n),AA(b,:)];
A_total(b,:)=[zeros(2,n),zeros(2,n)];
A_total(b+n,:)=[zeros(2,n),zeros(2,n)];
    case {'cl'}
b1=find(x==0);
S_total(b1,:)=[AA(b1,:),zeros(1,n)];
S_total(b1+n,:)=[zeros(1,n),AA(b1,:)];
A_total(b1,:)=[zeros(1,n),zeros(1,n)];
A_total(b1+n,:)=[zeros(1,n),zeros(1,n)];
b2=find(x==L);
S_total(b2,:)=[k*G*A*DAA(b2,:),k*G*A*AA(b2,:)];
S_total(b2+n,:)=[zeros(1,n),E*I*DAA(b2,:)];
A_total(b2,:)=[zeros(1,n),zeros(1,n)];
A_total(b2+n,:)=[zeros(1,n),zeros(1,n)];
end

%%



[lambda_vec,lambda]=eig(S_total,A_total);
[V,D]=eig(S_total,A_total,'qz');SS5=S_total*V-A_total*V*D;

T=S_total*lambda_vec-lambda*A_total*lambda_vec;




lambda=diag(lambda,0); 
[lambda,indice]=sort(lambda);
lambda_vec=lambda_vec(:,(indice(:)));

eigval = lambda(1,1);eigvec = lambda_vec(:,1);
S_total*eigvec - eigval*A_total*eigvec

% lambda=real(lambda);
%  lambda(lambda(:)==-inf | lambda(:)<=0)=NaN;
% 
% lambda=sqrt(lambda)*L^2*sqrt(rho*A/(E*I));%lambda=sqrt(lambda);




m=2;
sol_exacta=(m*pi/L)^2*sqrt((E*I)/(rho*A))*sqrt(1-(((m*pi/L)^2*E*I)/(k*G*A+(m*pi/L)^2*E*I)));
sol_exacta_norm=sol_exacta*L^2*sqrt(rho*A/(E*I));
%-------------graficos 


p=2;
lambda_mode_w(1:n,p)=lambda_vec(1:n,p)'*AA;
lambda_mode_phi_x(1:n,p)=lambda_vec(n+1:end,p)'*AA;
lambda_mode=[lambda_mode_w;lambda_mode_phi_x];


% SS3=lambda_vec(:,1)'*A_total*lambda_vec(:,1);
% % SS=*lambda_vec(:,p)-lambda(p,p)*A_total*lambda_vec(:,p);
% SS=(S_total+lambda(1)*A_total)*lambda_vec(:,1);
% 
% SS4=S_total*real(lambda_vec3)-A_total*real(lambda_vec3)*real(lambda3);
% SS4=S_total*(lambda_vec3)-A_total*(lambda3);
figure(1)
subplot(1,3,1);plot(x, lambda_mode_w(:,p));hold on;title(['w(' num2str(m) ')_{exact} = ' num2str(sol_exacta/(2*pi),'%6.4f')]);legend(['w(' num2str(m) ')=' num2str(sqrt(lambda(m))/(2*pi),'%6.4f')])
subplot(1,3,2);plot(x, lambda_mode_phi_x(:,p));hold on;title(['w(' num2str(m) ')_{exact} = ' num2str(sol_exacta_norm,'%6.4f')]);legend(['w(' num2str(m) ')=' num2str(lambda(m),'%6.4f')])

%%
freq=sqrt(lambda(m))/(2*pi);
vetor_carga=zeros(2*n,1); vetor_carga(2:n-1,1)=carga;
% sol_x=S_total\vetor_carga;
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sol_x=[lambda_mode_w(:,p);lambda_mode_phi_x(:,p)];
C=zeros(2*n,2*n);
% C=0.5*A_total+0.00001*S_total;
vetor_f=zeros(2*n,1);
%NEWMARK
% x_0=zeros(2*n,1);
x_0=sol_x;
v_0=zeros(2*n,1);
a_0=pinv(A_total)*(vetor_f-C*v_0-S_total*x_0); 
delta=1/2; alpha=1/3.9;
% delta_t=1/freq/100;   %delta t
delta_t=1/100000;
a0=1/(alpha*delta_t^2); a1=delta/(alpha*delta_t); a2=1/(alpha*delta_t); a3=1/(2*alpha)-1;
a4=delta/alpha-1; a5=(delta_t/2)*(delta/alpha-2); a6=delta_t*(1-delta); a7=delta*delta_t;

K_efe=S_total+a0*A_total+a1*C;

% t_final=10/freq;   % t final
t_final=0.1;
n_t=int16(t_final/delta_t+1);
t=zeros(n_t,1);
x_t=zeros(2*n,n_t); x_t(:,1)=x_0;
v_t=zeros(2*n,n_t); v_t(:,1)=v_0;
a_t=zeros(2*n,n_t); a_t(:,1)=a_0;
for i=2:n_t
  t(i)=t(i-1)+delta_t;
  F_efe=A_total*(a0*x_t(:,i-1)+a2*v_t(:,i-1)+a3*a_t(:,i-1))+C*(a1*x_t(:,i-1)+a4*v_t(:,i-1)+a5*a_t(:,i-1));
  x_t(:,i)=K_efe\F_efe;
  a_t(:,i)=a0*(x_t(:,i)-x_t(:,i-1))-a2*v_t(:,i-1)-a3*a_t(:,i-1);
  v_t(:,i)=v_t(:,i-1)+a6*a_t(:,i-1)+a7*a_t(:,i);
end

%interpolacao
x_w=zeros(n,n_t);
x_max=zeros(n_t,1);
for i=1:n_t
aux=x_t(1:n,i)'*AA;
x_w(:,i)=aux';
x_max(i)=x_w(ceil(n/4),i);   %!!!!!!!!!!!!!!!!!!!!!MUDAR DEPENDENDO DO MODO NATURAL!!!!!!!!!!
end

figure(2)
y_harm=x_max(1)*cos(2*pi*freq*t);
plot(t,y_harm);
hold on
plot(t,x_max);
hold on

ntt =double(2^nextpow2(n_t));
% X=fft(x_max,ntt);
X=fft(x_max);
X_mag=abs(X(1:ceil(n_t/2)));
[pk_vals, pk_locs]=findpeaks(X_mag);
% %remove peaks below threshold
% inds=find(X_mag(pk_locs)<1);
% pk_locs(inds)=[];

%determine frequencies
pk_freqs=zeros(length(pk_locs),1);
for i=1:length(pk_locs)
pk_freqs(i)=(pk_locs(i)-1)/t_final;
end

figure(3)
plot(X_mag);
hold on