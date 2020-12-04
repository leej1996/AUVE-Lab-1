N = 200;
nb = 20;

sigma1 = 2;
sigma2 = 5;
rho = 0.9;

% Get covariance matrix
Cxy = sigma1*sigma2*rho;
C = [sigma1^2, Cxy; Cxy, sigma2^2];

% Generate N realizations
x = chol(C, 'lower')*randn(length(C), N);

% Estimate mean and variance
m = mean(x,2)
cv = cov(x')

% Plot realizations
figure(1);
plot(x(1,:), x(2,:), 'x')

% Plot 91% confidence ellipses - real covariance (red) vs. estimated one (green)
hold on;
P0 = 0.91;
theta = linspace(0, 2*pi, 100);
X = sqrt(-2*log(1-P0))*chol(C, 'lower')*[cos(theta); sin(theta)];
plot(X(1,:), X(2,:), 'r');

hold on;
X = sqrt(-2*log(1-P0))*chol(cv, 'lower')*[cos(theta); sin(theta)] + m*ones(1, length(theta));
plot(X(1,:), X(2,:), 'g');