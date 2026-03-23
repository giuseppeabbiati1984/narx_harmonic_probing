close all
clear
clc

% this script shows how to compute:
% - the frequency combs of output and residual for an input frequency bin
% - the checks on the identifiability of the GRF frequency bin

% order of the GFRF
n = 2 ;

% input comb
wu = sym('wu',[n,1]) ;

%% output comb
wy = wu ;
if n > 2
    for i = 1:1:n-2
        wy = unique(kronsum(wy,wu),'rows') ;
    end
end
wy = [wy;sum(wu)] ;
wy = unique(wy,'rows') ;

%% residual comb
we = wy ;
if n > 1
    for i = 1:1:n-1
        we = unique(kronsum(we,wy),'rows') ;
    end
end
we = unique(we,'rows') ;

% max(abs(we))

%% check

iq = 0:1:n ;

cmd = 'combinations(iq';
for i = 1:1:n-1
    cmd = [cmd ',iq'] ;
end
cmd = [cmd,');'] ;

iq_vec = eval(cmd) ;
iq_vec = iq_vec.Variables ;

% retain only the rows that sum to m
iq_vec = iq_vec(sum(iq_vec,2)==n,:) ;

% remove the row equal to ones(1,m)
iq_vec(sum(iq_vec==ones(1,n),2)==2,:) = [] ;

% check to verify that the GRF bins are observable
check = iq_vec*wu==sum(wu) ;

%% codegen

vars = {wu} ;

matlabFunction(wy,'file',['eval_wy_' num2str(n)],'vars',vars) ;
matlabFunction(we,'file',['eval_we_' num2str(n)],'vars',vars) ;
matlabFunction(check,'file',['eval_check_' num2str(n)],'vars',vars) ;

return

%%

close all
clear
clc

dw = 0.1 ;
wb = 10 ;
w = -wb:dw:(wb-dw) ;

n = 4 ;

% select a randon frequency bin
wu = w(randi(numel(w),[n,1])).' ;

% if bins are distringuishable
if sum(eval_check(wu))==0
    % evaluate combs
    wy = eval(['eval_wy_' num2str(n) '(wu);']) ;
    we = eval(['eval_we_' num2str(n) '(wu);']) ;
end

%%

% test of the kronecker sum function
a = sym('a',[3,1]) ;
b = sym('b',[2,1]) ;
c = kronsum(a,b) ;

function c = kronsum(a,b)

Ia = ones(size(a)) ;
Ib = ones(size(b)) ;

c = kron(a,Ib)+kron(Ia,b) ;

end