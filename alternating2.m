% Load the matrices from .mat files
matrices = load('SystMatr.mat'); 
Mk = matrices.Mk;
A0 = matrices.A0;
B = matrices.B;
C = matrices.C;


% number of modules 
n = 2;


%% optimization 
k = 20000 * ones(1, n-1); % initialize guess for ki parameters 

k_nom = 10000 * ones(1, n-1);
Q_vals = {};
kivals = [];

J = 10 %number of refining iterations
y = 12000; % start with some bound for \gamma (H inf norm)
dy = 1; % push down guess for \gamma to try to get tighter bound
%%
ki = k;
% step 1: find feasible Q given values for \gamma, ki,
for jj = 1:10
ok = true; 
count = 0;
while ok && count<30
A = A0;
for i = 1:size(Mk, 1)
    % Multiply each layer by its corresponding scalar
    scaled_matrix = ki(i) * squeeze(Mk(i, :, :));
    
    % Add the scaled matrix to the sum matrix
    A = A + scaled_matrix;
end
disp('Matrix A:');
disp(A);
cvx_begin sdp
variable Q(2*n,2*n) semidefinite

X11 = (A')*Q + Q*A;
X = [X11, Q*B, C';...
B'*Q, -y*eye(n), zeros(n,n);...
C, zeros(n,n), -y*eye(n)];
minimize 1
subject to
-X == semidefinite(4*n,4*n)
cvx_end
if isequal(cvx_optval, Inf)
ok = false;
else %optimizing happens here
Q_last = Q;
y = y - dy; count = count + 1
end
end
disp('next step')
if any(any(isnan(Q)))
Q = Q_last; 
y = y + dy;
end
disp(y)
ok = true; 
% step 2: fix Q and update k1 
k1_last = ki; 
count = 0;
while ok && count < 200
cvx_begin sdp
variable ki(1,n-1) nonnegative
A = A0;
for i = 1:size(Mk, 1)
    % Multiply each layer by its corresponding scalar
    scaled_matrix = ki(i) * squeeze(Mk(i, :, :));
    
    % Add the scaled matrix to the sum matrix
    A = A + scaled_matrix;
end
X11 = (A')*Q + Q*A
X = [X11, Q*B, C';...
B'*Q, -y*eye(n), zeros(n,n);...
C, zeros(n,n), -y*eye(n)];
minimize 1%max(max(ki-k_nom))
subject to
-X == semidefinite(4*n,4*n)
cvx_end
if isequal(cvx_optval, Inf)
ok = false;
else
k1_last = ki; 
y = y - dy; count = count + 1
end
end
if any(isnan([ki]))
ki = k1_last;
y = y + dy;
end
disp(y)
disp(count)
end
%%
% verify true H_inf norm for parameters obtained from optimization
% procedure:

G = ss(A, B, C, zeros(n,n));
hnorm_calculated = hinfnorm(G);
