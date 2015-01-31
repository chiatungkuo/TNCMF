function out = couplenmf(S, P, D, R, beta, epsilon, maxitr, tol)
% solves the following constrained optimization
%   minimize ||S - Lp * C * Ld'||^2 + beta * (||P - Lp * Tp'||^2 + ||D - Ld * Td'||^2) 
%   subject to 
%       0 <= Lp, Ld <= 1, 0 <= Tp, Td, 0 <= C <= max(S),
%       Lp(k, i) * Lp(k, j) <= epsilon for all k, for i != j
%       and Ld(k, i) * Ld(k, j) <= epsilon for all k, for i != j
%
% Inputs
%   S: n by n matrix
%   P: n by t matrix
%   D: n by t matrix
%   R: desired dimension of the output core 
%   beta: tradeoff paramenter in the objective function (default: 1)
%   epsilon: constant between 0 and 1 specifying the diversity constraint
%       use [] for no constraints (default: [])
%   maxitr: maximum of iterations allowed (default: 100)
%   tol: convergence criterion where the algorithm stops when the relative
%       change in objective is smaller than tol between two iterations
%       (default: 1e-4)
%
% Outputs
%   out: cell array such that out{1} = Lp, out{2} = Ld, out{3} = C,
%       out{4} = Tp and out{5} = Td
%

if ~exist('beta', 'var'), beta = 1; end
if ~exist('epsilon', 'var'), epsilon = []; end
if ~exist('maxitr', 'var'), maxitr = 100; end
if ~exist('tol', 'var'), tol = 1e-4; end

normS = norm(S, 'fro');
normP = norm(P, 'fro');
normD = norm(D, 'fro');

maxS = max(S(:));

% Initialize sparse random factor and core matrices
density = 0.2;
Lp = sprand(size(S, 1), R(1), density);
Ld = sprand(size(S, 2), R(2), density);
coreS = rand(R(1), R(2));

if isempty(P), Tp = []; else Tp = sprand(size(P, 2), R(1), density); end
if isempty(D), Td = []; else Td = sprand(size(D, 2), R(2), density); end

Lp = full(Lp);
Ld = full(Ld);

out = cell(1, 5);

objval = normS^2 + beta*(normP^2 + normD^2);

% Main loop
for iter = 1:maxitr
    
    objold = objval;
    
    % Update core constraining entries between 0 and 1 (inclusive)
    A = kron(Ld, Lp);
    b = reshape(S, [numel(S) 1]);
    
    if ~isempty(epsilon) && epsilon > 0
        % solve using coordinate descent with lower/upper bounds
        coreStmp = reshape(coreS, [R(1)*R(2) 1]);  
        coreStmp = coordupdate(A', b', coreStmp', 0, maxS, []);    
        coreS = reshape(coreStmp(:), R);
    else
        % solve using NNLS by Kim & Park
        [coreStmp, ~, ~, succ] = nnlsm_blockpivot(A'*A, A'*b, 1);
        if succ == 1, coreS = reshape(coreStmp, R); else fprintf('NNLS subproblem ill-conditioned.\n'); continue; end
    
        % solve using NNLS by Bro
%       coreS = reshape(fastnnls(A'*A, A'*b), R);
    end
       
    % Update Lp     
    A = [coreS*Ld' sqrt(beta)*Tp'];
    B = [S sqrt(beta)*P];
    
    if ~isempty(epsilon) && epsilon > 0
        Lp = coordupdate(A, B, Lp, 0, 1, epsilon);
    else
        [Lptmp, ~, ~, succ] = nnlsm_blockpivot(A*A', A*B', 1);
        if succ == 1, Lp = Lptmp'; else fprintf('NNLS subproblem ill-conditioned.\n'); continue; end

        % normalize Lp and put weights into core
        lambda = sqrt(sum(Lp.^2))';
        Lp = sparse(Lp)*spdiags(1./lambda, 0, R(1), R(1));
        coreS = spdiags(lambda, 0, R(1), R(1))*coreS;
    end
                
    % Update Ld
    A = [coreS'*Lp' sqrt(beta)*Td'];
    B = [S' sqrt(beta)*D];
    
    if ~isempty(epsilon) && epsilon > 0
        Ld = coordupdate(A, B, Ld, 0, 1, epsilon);
    else
        [Ldtmp, ~, ~, succ] = nnlsm_blockpivot(A*A', A*B', 1);
        if succ == 1, Ld = Ldtmp'; else fprintf('NNLS subproblem ill-conditioned.\n'); continue; end
    
        % normalize Ld and put weights into core
        lambda = sqrt(sum(Ld.^2))';
        Ld = sparse(Ld)*spdiags(1./lambda, 0, R(2), R(2));
        coreS = coreS*spdiags(lambda, 0, R(2), R(2));  
    end

    % Update Tp
    if ~isempty(Tp) 
        [Tptmp, ~, ~, succ] = nnlsm_blockpivot(Lp'*Lp, Lp'*P, 1);
        if succ == 1, Tp = Tptmp'; else fprintf('NNLS subproblem ill-conditioned.\n'); continue; end
    end
    
    % Update Td
    if ~isempty(Td)
        [Tdtmp, ~, ~, succ] = nnlsm_blockpivot(Ld'*Ld, Ld'*D, 1);
        if succ == 1, Td = Tdtmp'; else fprintf('NNLS subproblem ill-conditioned.\n'); continue; end
    end
    
    % check change in objective value
    out{1} = Lp; out{2} = Ld; out{3} = coreS; out{4} = Tp; out{5} = Td;
    
    Shat = Lp*coreS*Ld';
    Sfit = 1 - norm(S - Shat, 'fro')/normS;    
    if ~isempty(P)
        Phat = Lp*Tp';
        Pfit = 1 - norm(P-Phat, 'fro')/normP;
    else
        Phat = []; 
        Pfit = 0;
    end
    if ~isempty(D)
        Dhat = Ld*Td';
        Dfit = 1 - norm(D-Dhat, 'fro')/normD;
    else
        Dhat = [];
        Dfit = 0;
    end  
    
    objval = norm(S-Shat, 'fro')^2 + beta*(norm(P-Phat, 'fro')^2 + norm(D-Dhat, 'fro')^2);
    objchange = (objold - objval)/objold;
    
    fprintf('Iter %2d: obj. = %6.4e\tSfit = %6.4e\tPfit = %6.4e\tDfit = %6.4e\tobj. delta = %6.4e\n', iter, objval, Sfit, Pfit, Dfit, objchange);
    
    % Check for convergence
    if (iter > 1) && (objchange < tol)
        break;
    end
end

end



function X = coordupdate(A, B, X0, lb, ub, epsilon)
% Solves the constrained problem
%   minimize ||B - XA||_F^2
%   subject to lb <= X <= ub
% entry by entry while fixing all other constant (coordinate descent)

dim = size(A, 1);
for k = 1:dim
    if k == 1
        res = B - X0(:, 2:end)*A(2:end, :);
    else
        res = res - sol*A(k-1, :) + X0(:, k)*A(k, :);
    end

    x = A(k, :)';
    denom = x'*x;
    if denom < eps, continue; end
    
    sol = res*x/denom;
    sol = max(sol, lb);
    sol = min(sol, ub);
    if ~isempty(epsilon) && epsilon > 0
        diversityUpper = min(epsilon./X0(:, [1:k-1, k+1:end])')';
        sol = min(sol, diversityUpper);
    end

    X0(:, k) = sol;
end
X = X0;
end
