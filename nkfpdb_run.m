%MODEL_RUN.M
%Luisa Lambertini
%Date November 2013
close all
clc
clear all
global s2 t2 approx
%degree of approximation in the linearization around the steady state
approx=1;
[fx,fxp,fy,fyp,f] = nkfpdb_model;

%Steady State and Parameter Values
%original: [gammma,rhoz_star,betta,siggma,varphi,allpha,nu,theta,epsilon,phipi,phiy,rho,delta_m,pssi,int,infl,inflf,b,bl,bf,bss,mcss,mcf,M,Mss,lambda,lambdaf,z_star,P_z,gov,govss] = nkfpdb_ss;
[gammma,rhoz_star,betta,siggma,varphi,allpha,nu,theta,epsilon,phipi,phiy,rho,delta_m,pssi,int,infl,inflf,...
    b,bl,bf,bss,mcss,mcf,M,Mss,lambda,lambdaf,z_star,P_z,gov,govss,c,cf,profit, profitf, tauf, coef1, coef2, coef3,...
    xi, Iota, Iotaf,Mf,z,w,n,tau,nnf,yy,yyss,mc,z_starf,intf,zf,wf,yyf,Ml] = nkfpdb_ss;

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx_new(nfy,nfx,nfyp,nfxp,1);
x=[z_star, M];
y=[int, infl, c, z, lambda, yy, n, w, mc, profit, tau, Iota];
tt=sort(abs(diag(s2).\diag(t2)));
nx=size(x,2);
tt(nx)
x0 =[20, 0];

IR=ir(gx,hx,x0,40);
IR(:,1)=exp(int)*4*(1+IR(:,1)/100);
IR(:,2)=4*IR(:,2);

% %% add more graph titles as needed an in the same order as your variables
% Leave Interest rate and inflation as the first 2 title graphs
% %Make sure to have Interest rate and inflation as your first 2 variables
% in your "y" vector
titlegraph=['Int Rate    ';...
            'Inflation   ';...
            'Consumption ';...
            'Emission(z) ';...
            'Lambda      ';...
            'Output(y)   ';...
            'Labor(n)    ';...
            'Real wage(w)';...
            'mc          ';...
            'Profit      ';...
            'Tax(tau)    ';...
            'Iota        ';...
            'World Z*    ';...
            'Pollution(M)'];

nrow = 4;  %number of rows in the plot;
ncol = 4;  %number of column in the plot;
for j=1:size(IR,2)
    subplot(nrow,ncol,j);
    plot(IR(:,j))
    title(titlegraph(j,:));
      xlabel('Quarters after 0');
      if j==1
          ylabel('level, ann'); 
      elseif j==2
          ylabel('%dev, ann');
      else
      ylabel('% dev from stst'); 
      end
end  