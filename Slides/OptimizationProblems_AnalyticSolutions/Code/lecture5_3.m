clc; clear; close all;

v = randn(2, 1); v = v / norm(v);
c = randn(2, 1);
a = randn(2, 1);


M = [eye(2),      zeros(2, 1),-eye(2);
     zeros(1, 2), zeros(1, 1), v';
    -eye(2),      v,           zeros(2, 2)];
vec = [a; 0; -c];

rank(M)
size(M)

var = pinv(M) * vec;

x_analytic = var(1:2);

%%%%%%%%%%%%

ksi2 = v'*(a-c);
x_analytic2 = v*v'*a + (eye(2) - v*v')*c;
lambda2 = (eye(2) - v*v')*(c - a);
var2 = [x_analytic2; ksi2; lambda2];

disp("norm(var), norm(var2), norm(var-var2)")
[norm(var), norm(var2), norm(var-var2)]

%%%%%%%%%%%%
% check

cvx_begin
variables x(2);
variables al(1);

minimize( 0.5*(x - a)'*(x - a) );
subject to
    x == v*al + c;
cvx_end

x_cvx = x;


%%%%%%%%%%%%%%%%%%%%%%%%

w = null(v');

var2 = pinv([v, w]) * (a - c);

alpha = var2(1);


x_geom = v*alpha + c;


%%%%%%%%%%%%%%%%%%%%%%%%



disp(" ")
disp(" ")
disp("x_analytic, x_cvx, x_analytic-x_cvx, x_geom-x_cvx")
[x_analytic, x_cvx, x_analytic-x_cvx, x_geom-x_cvx]

disp("norm(x_analytic), norm(x_cvx), norm(x_analytic-x_cvx), norm(x_geom-x_cvx)")
[norm(x_analytic), norm(x_cvx), norm(x_analytic-x_cvx), norm(x_geom-x_cvx)]


Count = 50;
res.X = zeros(Count, 1);
res.Y = zeros(Count, 1);

for i = 1:Count
    phi = (i-1)*2*pi/(Count-1);
    res.X(i) = cos(phi);
    res.Y(i) = sin(phi);
end



scale = 3;

figure;
plot([c(1) - scale*v(1); c(1) + scale*v(1)], [c(2) - scale*v(2); c(2) + scale*v(2)], '-', ...
    'LineWidth', 3, 'Color', 'k'); hold on
plot(a(1), a(2), 'o', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r'); hold on

plot(x_cvx(1), x_cvx(2), 's', 'LineWidth', 0.5, 'MarkerSize', 15, 'MarkerFaceColor', 'b'); hold on
plot([x_cvx(1); a(1)], [x_cvx(2); a(2)], '-', 'LineWidth', 2, 'Color', 'b'); hold on



plot(x_analytic(1), x_analytic(2), 'o', 'LineWidth', 1, 'MarkerSize', 20, 'Color', 'b'); hold on
plot([x_analytic(1); a(1)], [x_analytic(2); a(2)], '--', 'LineWidth', 4, 'Color', 'r'); hold on


plot(x_geom(1), x_geom(2), 'o', 'LineWidth', 0.5, 'MarkerSize', 30, 'Color', 'k'); hold on



plot(0.4*res.X + a(1), 0.4*res.Y + a(2), '-', 'LineWidth', 0.5, 'Color', [0.2, 0.2, 0.2]); hold on
plot(res.X + a(1),         res.Y + a(2), '-', 'LineWidth', 0.5, 'Color', [0.2, 0.2, 0.2]); hold on
plot(1.6*res.X + a(1),     1.6*res.Y + a(2), '-', 'LineWidth', 0.5, 'Color', [0.2, 0.2, 0.2]); hold on


axis equal;