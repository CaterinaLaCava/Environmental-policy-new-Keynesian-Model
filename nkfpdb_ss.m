function [gammma,rhoz_star,betta,siggma,varphi,allpha,nu,theta,epsilon,phipi,phiy,rho,delta_m,pssi,int,infl,inflf,...
    b,bl,bf,bss,mcss,mcf,M,Mss,lambda,lambdaf,z_star,P_z,gov,govss,c,cf,profit, profitf, tauf, coef1, coef2, coef3, xi, Iota, Iotaf, Mf,z,w,n,tau,nnf,yy,yyss,mc,z_starf,intf,zf,wf,yyf,Ml] = nkfpdb_ss
%added c,cf,profit, profitf, tauf,coef1, coef2, coef3, xi, Iota,
%Iotaf,Mf,z,w,n,tau,nnf,yy,yyss,mc,z_starf,intf,zf,wf,yyf
global siggma allpha varphi epsilon nu betta bss z P_z pssi delta_m M rhoz_star 
global phipi phiy gammma rho bigtheta rho theta int infl b mc coef1 coef2 coef3 xi z_star %added
%%%% Check for any missing parameter!
betta = 0.99;
siggma=1;
varphi = 1;
allpha = 0.3;
nu = 1;
theta = 2/3;
epsilon = 11;
phipi = 1.5;
phiy = 0.5/4;
rho = 0.9;
rhoz_star = 0.8;
delta_m = 1 - 0.9979;
pssi = 0.4;                   %%% you are asked to find pssi%added: Bahrain; 
bigtheta=(1-allpha)/(1-allpha+allpha*epsilon);
gammma=(1-theta*betta)*(1-theta)*bigtheta/theta;
coef1 = 1/1000;                                                 %added
coef2 = 6/100;                                                  %added
coef3 = 6/100000;                                               %added
xi = 0.1; %added

%added comment: % steady state values: non constant so shouldnt be global??
P_z = 0.5;
int = 1/betta - 1;
infl=1;
b=1; 
M = 41.912; %%% you are asked to find M, divide its value by 10 %added 419.12/10
Mss = M;
z_star = 1;
mc=(epsilon-1)/epsilon;
mcss=mc;

%%% Steady state with fsolve (generic solution)
guess=[0.27 0.19 0.34 1 1 0.2 0.2 0.8 0.6];
options0=optimset('MaxFunEvals',300000,'MaxIter',500000,'TolFun',1e-12);
xxx=fsolve('steady_state_nkfpdb',guess,options0);
steady_state_nkfpdb(xxx)

%%%% Extract the steady state variables from the "xxx" vector
%copied same order we have in y = [...] % variables already defined in this file:infl mc int 
c=xxx(1);
z=xxx(2);
lambda=xxx(3);
yy=xxx(4);
n=xxx(5);
w=xxx(6); 
profit=xxx(7);
tau=xxx(8);
Iota=xxx(9); 


%%%%%

yyss=yy;
gov = 0.4*yyss;
govss = gov;
bss = b;

%take logs
infl = log(infl);
mc = log(mc);
lambda = log(lambda);
b = log(b);
c=log(c); %added
z=log(z); %added
n=log(n); %added
yy=log(yy); %added
w=log(w); % added
int=log(int); %added
profit=log(profit);%added
tau=log(tau);%added
Iota=log(Iota); %added
z_star=log(z_star);%added
M=log(M);%added
Mss = log(Mss); %added


%define t+1 and t-1 variables at the steady state
lambdaf=lambda;
mcf=mc;
inflf=infl;
bl = b;
bf = b;
Mf = M; %added
Ml = M; %added
z_starf = z_star; %added
tauf=tau; %added
cf = c; %added
profitf=profit; %added
Iotaf=Iota; %added
nnf=n;%added
yyf=yy;%added
wf=w;%added
zf=z;%added
intf=int;%added


