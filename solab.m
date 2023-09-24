%
% Function: solab
%
% Purpose: Solves for the recursive representation of the stable solution 
% to a system of linear difference equations.
%
% Inputs: Two square matrices a and b and a natural number nk
%
% a and b are the coefficient matrices of the difference equation
%
% a*x(t+1) = b*x(t)
% 
% where x(t) is arranged so that the state variables come first, and
%
% nk is the number of state variables.
%
% Outputs: the decision rule f and the law of motion p. If we write
%
% x(t) = [k(t);u(t)] where k(t) contains precisely the state variables, then
% 
% u(t)   = f*k(t) and
%
% k(t+1) = p*k(t).
%
% Calls: reorder
%

function [f,p] = solab(aa,b,nk);
global sv t z11 z21 s2 t2 nk
[sv,t,q,zz] = qz(aa,b);  % upper triangular factorization of the matrix pencil b-za
s2=sv;
t2=t;
[sv,t,q,zz] = reorder(sv,t,q,zz);   % reordering of generalized eigenvalues in ascending order

z21 = zz(nk+1:end,1:nk);
z11 = zz(1:nk,1:nk);

if rank(z11)<nk;
 error('Invertibility condition violated')
end

z11i = z11\eye(nk);
s11 = sv(1:nk,1:nk);
t11 = t(1:nk,1:nk);

if abs(t(nk,nk))>abs(sv(nk,nk)) | abs(t(nk+1,nk+1))<abs(sv(nk+1,nk+1));
 warning('Wrong number of stable eigenvalues.');
end

dyn = s11\t11; %This is inv(s11)*t11

f = real(z21*z11i);      % The real function takes away very small imaginary parts of the solution
p = real(z11*dyn*z11i);