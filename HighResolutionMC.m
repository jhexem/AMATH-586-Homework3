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
    U(2, k) = U(2, k-1) - (dt / (2 * dx)) * (U(3, k-1) - U(1, k-1)) + ...
            (dt^2 / (2 * dx^2)) * (U(3, k-1) - 2 * U(2, k-1) + U(1, k-1));
    for j = 3:(length(xvals)-1)
        theta1 = (U(j, k-1) - U(j-1, k-1)) / (U(j+1, k-1) - U(j, k-1));
        theta2 = (U(j-1, k-1) - U(j-2, k-1)) / (U(j, k-1) - U(j-1, k-1));
        theta1vec = [(1+theta1)/2, 2, 2*theta1];
        theta2vec = [(1+theta2)/2, 2, 2*theta2];
        F1 = U(j, k-1) + (1/2) * (1-(dt/dx)) * (U(j+1, k-1) - U(j, k-1)) * ...
            max(0, min(theta1vec));
        F2 = U(j-1, k-1) + (1/2) * (1-(dt/dx)) * (U(j, k-1) - U(j-1, k-1)) * ...
            max(0, min(theta2vec));
        U(j, k) = U(j, k-1) - (dt / dx) * (F1 - F2);
    end
    U(end, k) = (1 - dt/dx) * U(end, k-1) + (dt/dx) * U(end-1, k-1);
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
title("High Resolution Scheme with MC Flux Limiter Solution at t=0.8")
legend({"True Solution", "High Resolution Scheme Solution"}, 'Location', 'northwest')
xlabel("x")
ylabel("u(0.8, x)")