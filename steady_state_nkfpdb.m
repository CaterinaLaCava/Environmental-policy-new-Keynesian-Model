function eq = steady_state_nkfpdb(xxx) 
global siggma allpha varphi epsilon nu P_z pssi delta_m M 
global z_star           %added
%%%%% Define the variables in vector "xxx" to extract the steady 
%%%%% state variables after solving with fsolve
%%%%% should be the same as what you define in the "_ss" file
c=xxx(1);
z=xxx(2);
lambda=xxx(3);
yy=xxx(4);
n=xxx(5);
w=xxx(6); %already defined in this file:infl mc int 
profit=xxx(7);
tau=xxx(8);
Iota=xxx(9); 

%%%%
%y=[c z lambda yy n w profit tau Iota]

eq(1) = 1/c^(1/siggma) - lambda; %foc C at steady state                         %added
eq(2) = -yy+c+0.5*yy; % Market clearing
eq(3) = -yy+(1-Iota)*n^(1-allpha); %% production function                       %added
eq(4) = -w+(epsilon-1)/epsilon*(1-allpha)*yy/n; %% Marg cost
eq(5) =  - n^(1/varphi)*nu - lambda*w*(tau - 1); %% FOC wrt labor                %added
eq(6) = -0.5*yy+tau*w*n  +P_z*z; 
eq(7) = -z+pssi*yy; %% define Z                                                 %added
eq(8) = -M + (1-delta_m)*M + z ;%+ 0.1*z_star;  %% law of motion of M at steady state%added
eq(9) = -profit+yy-w*n-P_z*z; %%% define profit                                  %added
