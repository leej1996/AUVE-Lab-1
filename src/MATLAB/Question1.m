%% Main function for the questions of the first part
function [] = Question1()

% Solves Question 1.1
Question11();

% Solves Questions 1.2
Question12();

end

%% Solves Question 1.1
function [] = Question11()

% Parameters
m = 0.95; % mean
p = 0.5; % prob
sigma = sqrt(1 - (m^2)); % sigma

% For N=100 and N=4000 realizations
for N = [100 4000]
    
    %% X normally distributed
    Xn= randn(N,1);
    Yn=pdf('Normal',-5:0.01:5,0,1);
    figure;
    histogram(Xn,'Normalization','pdf');
    hold on;
    plot(-5:0.01:5,Yn);
    title(['N = ',num2str(N),' Normal']);
    meann=mean(Xn);
    stdn=std(Xn);
    disp('***** Normal *****');
    disp(['Mean: ',num2str(meann),' Std: ',num2str(stdn)]);
    
    figure;
    plot(1:N,Xn,'*');
    title(['Normal Independent Random Signal, N = ',num2str(N)]);
    
    %% X uniformly ditributed
    Xu=2*sqrt(3)*(rand(N,1)-0.5);
    Yu=pdf('Uniform',-2:0.01:2,-sqrt(3),sqrt(3));
    figure;
    histogram(Xu,'Normalization','pdf');
    hold on;
    plot(-2:0.01:2,Yu);
    title(['N = ',num2str(N),' Uniform']);
    meanu=mean(Xu);
    stdu=std(Xu);
    disp('***** Uniform *****');
    disp(['Mean: ',num2str(meanu),' Std: ',num2str(stdu)]);
    
    figure;
    plot(1:N,Xu,'*');
    title(['Uniform Independent Random Signal, N = ',num2str(N)]);
    
    %% X driven by a Gaussian mixture
    Xs = randn(N,1)*sqrt(1-m*m)+m;
    k = find(rand(N,1)>0.5);
    Xs(k) = Xs(k) -2*m;
    Ys= p*pdf('Normal',-5:0.01:5,m,sigma)+(1-p)*pdf('Normal',-5:0.01:5,-m,sigma);
    figure;
    histogram(Xs,'Normalization','pdf');
    hold on;
    plot(-5:0.01:5,Ys);
    title(['N = ',num2str(N),' Symmetric']);
    means=mean(Xu);
    stds=std(Xu);
    disp('***** Symmetric *****');
    disp(['Mean: ',num2str(means),' Std: ',num2str(stds)]);
    
    figure;
    plot(1:N,Xs,'*');
    title(['Symmetric Independent Random Signal, N = ',num2str(N)]);
    
end

end

%% Solves Question 1.2
function [] = Question12()

% Parameters
sigma1 = 2; % Sigma X1
sigma2 = 5; % Sigma X2
p = 0.9; % correlation coefficient

C = [sigma1^2, sigma1*sigma2*p;sigma1*sigma2*p, sigma2^2]; % Covariance matrix
Xj = chol(C,'lower')*randn(2,200);

meanj= mean(Xj,2);
covj= cov(Xj');

P0 = 0.91;
theta = linspace(0,2*pi,100);
X = sqrt(-2*log(1-P0))*chol(C,'lower')*[cos(theta);sin(theta)] + meanj*ones(1,length(theta));

disp('***** Joint *****');
disp(['Mean: ',num2str(meanj')]);
disp('Covariance:');
disp(covj);

figure;
plot(Xj(1,:), Xj(2,:), '*');
hold on;
plot(X(1,:),X(2,:));
title('Joint Distribution with 91% Confidence Ellipse');

end