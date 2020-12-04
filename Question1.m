%% Question 1.1
m = 0.95;
nb = 16;
p = 0.5;
sigma = sqrt(1 - (m^2));

for N = [100 4000]
    X_n = randn(N,1);
    X_u = 2*sqrt(3)*(rand(N,1)-0.5);
    X_s = randn(N,1)*sqrt(1-m*m)+m;
    k = find(rand(N,1)>0.5);
    X_s(k) = X_s(k) -2*m;
    
    %Plot normal distribution
    [n,b] = hist(X_n,nb);
    figure;
    hold on;
    bar(b,n/(b(2)-b(1))/sum(n),1);
    title(['N = ',num2str(N),' Normal']);
    x = -5:0.01:5;
    pdf_normal = 1/sqrt(2*pi)*exp(-0.5*x.^2);
    plot(x,pdf_normal);
    mean_normal = mean(X_n);
    std_normal = std(X_n);
    fprintf('The mean and standard deviation of the Normal R.V. at N = %d is: %f, %f\n',N, mean_normal, std_normal);
    
    figure;
    plot(1:N,X_n,'*');
    title(['Normal Independent Random Signal, N = ',num2str(N)]);
    
    %Plot uniform distribution
    [n,b] = hist(X_u,nb);
    figure;
    hold on;
    bar(b,n/(b(2)-b(1))/sum(n),1);
    title(['N = ',num2str(N),' Uniform']);
    x = -2:0.01:2;
    pdf_uniform = (1/(2*sqrt(3)))*(x >= -sqrt(3) & x <= sqrt(3));
    plot(x,pdf_uniform);
    mean_uniform = mean(X_u);
    std_uniform = std(X_u);
    fprintf('The mean and standard deviation of the Uniform R.V. at N = %d is: %f, %f\n',N, mean_uniform, std_uniform);
    
    figure;
    plot(1:N,X_u,'*');
    title(['Uniform Independent Random Signal, N = ',num2str(N)]);
        
    %Plot Symmetric Gaussian Mixture
    [n,b] = hist(X_s,nb);
    figure;
    hold on;
    bar(b,n/(b(2)-b(1))/sum(n),1);
    title(['N = ',num2str(N),' Symmetric']);
    x = -5:0.01:5;
    pdf1 = p*(1/sqrt(2*pi)/sigma*exp(-0.5*((x-m)/sigma).^2));
    pdf2 = (1-p)*(1/sqrt(2*pi)/sigma*exp(-0.5*((x+m)/sigma).^2));
    pdf_symmetric = pdf1 + pdf2;
    plot(x,pdf_symmetric);
    mean_symmetric = mean(X_s);
    std_symmetric = std(X_s);
    fprintf('The mean and standard deviation of the Uniform R.V. at N = %d is: %f, %f\n',N, mean_symmetric, std_symmetric);
    
    figure;
    plot(1:N,X_s,'*');
    title(['Symmetric Independent Random Signal, N = ',num2str(N)]);
end

%% Question 1.2
sigma1 = 2;
sigma2 = 5;
p = 0.9;

C = [sigma1^2, sigma1*sigma2*p;sigma1*sigma2*p, sigma2^2];

x = chol(C,'lower')*randn(2,200);

m = mean(x,2);
covariance = cov(x');
disp("The mean is: ");
disp(m);
disp("The covariance is: ");
disp(covariance);

figure;
hold on;
plot(x(1,:), x(2,:), '*');

P0 = 0.91;
theta = linspace(0,2*pi,100);
X = sqrt(-2*log(1-P0))*chol(C,'lower')*[cos(theta);sin(theta)] + m*ones(1,length(theta));

plot(X(1,:),X(2,:));
title('Joint Distribution with 91% Confidence Ellipse');