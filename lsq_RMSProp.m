function [ X, resnorm, exitflag ] = lsq_RMSProp( f, X0, options)
%RMSPROP Summary of this function goes here
%   Detailed explanation goes here

% set constant
alpha = 1;%global descent rate
nita = 0.5; % momentum coefficient
gamma = 0.9; % decay rate

% Initialize
beta = randn(size(X0)).*X0;
v = randn(size(X0)).*X0;
r = 0;
if isfield(options, 'Display')
    switch options.Display
        case 'iter'
            IterDisp = true;
            FinalDisp = true;
        case 'final'
            IterDisp = false;
            FinalDisp = true;
        case 'none'
            IterDisp = false;
            FinalDisp = false;
    end
else % default as final
    IterDisp = false;
    FinalDisp = true;
end

if isfield(options, 'MaxIterations')
    MaxIter = options.MaxIterations;
else
    MaxIter = inf;
end

if isfield(options,'FunctionTolerance')
    res_tol = options.FunctionTolerance;
else
    res_tol = 1e-6;
end

if isfield(options,'StepTolerance')
    x_tol = options.StepTolerance;
else
    x_tol = 1e-6;
end

if IterDisp
    fprintf('                                         Norm of      First-order\n')
    fprtinf(' Iteration  Func-count     f(x)          step         optimality\n')
end

x = X0;
[y, dydx] = f(x);
g = 2*y'*dydx;

res = norm(y);
iter = 0;
fcount = 1;
if IterDisp
    fprintf(' %-9d  %-9d %-12.5g %-12.5g %-12.5g',iter,fcount,res,norm(v),norm(g));
end

while res>res_tol || norm(v)>x_tol || iter<=MaxIter
    beta = x+nita*v;
    [y, dydx] = f(beta);
    g = 2*dydx'*y;
    r = gamma*r+(1-gamma)*g'*g;
    v = nita*v - alpha*g/sqrt(r);
    x = x + v;
    res = norm(y);
if IterDisp
    fprintf(' %-9d  %-9d %-12.5g %-12.5g %-12.5g',iter,fcount,res,norm(v),norm(g));
end
end

X = x;
resnorm = res;
exitflag = 1;
end

