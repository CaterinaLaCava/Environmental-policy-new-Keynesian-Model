function [fx,fxp,fy,fyp,f] = nkfpdb_model
global approx
approx=1;
%Define parameters
syms betta siggma allpha delta_m varphi nu phipi phiy pssi epsilon rhoz_star tau tauf gammma P_z
syms gov theta coef1 coef2 coef3 xi govf %added
%Define variables 
syms c cf z zf z_star z_starf b bl bf bss lambda lambdaf yy yyf n nnf w wf 
syms infl inflf yyss mc mcf mcss int intf M Mf Mss Ml
syms Iota profit profitf Iotaf                                               %added

%define utility function
utility=c^(1-1/siggma)/(1-1/siggma)-nu/(1+1/varphi)*(n)^(1+1/varphi);
%define budget constraint of consumer
bc=-c-b/(1+int)+bl/infl+w*n*(1-tau)+profit;
bcf=-cf-bf/(1+intf)+b/inflf+wf*nnf*(1-tauf)+profitf;
%define lagrangean
lagrangian_cons=utility+lambda*bc+betta*lambdaf*bcf;
f1= jacobian(lagrangian_cons,c);                                              %added:%%% FOC wrt c
f2= jacobian(lagrangian_cons,b);                                              %added:%%% FOC wrt b
f3= jacobian(lagrangian_cons,n);                                              %added:%%% FOC wrt n
f4= jacobian(lagrangian_cons,lambda);                                         %added:%%% FOC wrt lambda
%define output
f5= -Iota+coef1+coef2*M+coef3*M^2;                                           %added:%%%Iota polynomial
f6= -yy+(1-Iota)*n^(1-allpha);                                               %added:%%%define output
f7= -z+pssi*yy;                                                              %added:%%% define z *Q then why are we using another symbolic variable z?
f8= -Mf + (1-delta_m)*M + z + 0.1*z_star;                                    %*Q: why is this 1 period in the future? A: matlab technicality
%governemnt budget constraint
f9= -gov +P_z*z + tau*w*n;%-gov-bl/infl+tau*w*n+b/(1+int)+P_z*z;             %added:%define Government BC
%monetary policy
f10= -int-log(betta) + phipi*log(infl)+ phiy*log(yy/yyss); 
%phillips curve
f11= -log(infl)+betta*log(inflf)+gammma*(log(mc) -log(mcss));
%Marginal Cost
f12= -mc+w/(((1-allpha)*n^(-allpha)*(1-Iota))); 
%exogenous shock to technology
f13= log(z_starf)-log(z_star)*rhoz_star; 
f14= -profit+yy-w*n-P_z*z;                                                   %added:%%% define profit
%create f
f= [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13;f14];
% Define the vector of controls, y, and states, x
x=[z_star, M];
y=[int, infl, c, z, lambda, yy, n, w, mc, profit, tau, Iota];                 %added
xp=[z_starf, Mf];
yp=[intf, inflf, cf, zf, lambdaf, yyf, nnf, wf, mcf, profitf, tauf, Iotaf];

lv=[x,y,xp,yp];
for i=1:size(lv,2)
f = subs(f,lv(i),exp(lv(i)));
end

%Compute analytical derivatives of f
[fx,fxp,fy,fyp]=anal_deriv(f,x,y,xp,yp,approx); 

