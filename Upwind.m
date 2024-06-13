u0 = @(x) (x >= 0 & x<=0.6) .* exp(-100 .* (x-0.3).^2) + (x > 0.6 & x <= 1);
uexact = @(x, t) u0(x - t);

N = 200;
dx = 1 / (N + 1);
dt = 0.8 * dx;

xvals = 0:dx:2;
tvals = 0:dt:1;

[T, X] = meshgrid(tvals, xvals);

sol1 = zeros(length(xvals), 1);  %solution at t=0.5
sol2 = zeros(length(xvals), 1);  %solution at t=0.8

U = zeros(length(xvals), length(tvals));
U(:, 1) = u0(xvals);

for k = 2:length(tvals)
    for j = 2:(length(xvals))
        U(j, k) = (1 - dt/dx) * U(j, k-1) + (dt/dx) * U(j-1, k-1);
    end
    if (tvals(k) >= 0.5 && norm(sol1) == 0)
        sol1 = U(:, k);
    end
    if (tvals(k) >= 0.8 && norm(sol2) == 0)
        sol2 = U(:, k);
    end
end


plot(xvals, uexact(xvals, 0.8), "LineWidth", 2)
hold on;
plot(xvals, sol2, "LineWidth", 2)
title("Upwind Scheme Solution at t=0.8")
legend({"True Solution", "Upwind Scheme Solution"}, 'Location', 'northwest')
xlabel("x")
ylabel("u(0.8, x)")
