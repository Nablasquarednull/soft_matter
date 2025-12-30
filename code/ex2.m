clear all; close all; clc;
% ANALYSIS OF SIMULATION OF 3-D FJC MODEL OF POLYMERS

% PARAMETERS
b = 3.0 ;
N = 100 ;
T = 1000 ;
bins = 50;

% INITIALIZATION
X = zeros(T,N+1) ;
Y = zeros(T,N+1) ;
Z = zeros(T,N+1) ;

% READ INPUT FILE
filename = sprintf('simulation_FJC_b=%.1f_N=%d_T=%d.xyz', b, N, T);
fid = fopen(filename, 'r');

for t = 1:T
    fgetl(fid) ;
    fgetl(fid) ;

    for i = 1:N+1
        line = fgetl(fid) ;
        data = strsplit(line) ;
        X(t,i) = str2double(data(2)) ;
        Y(t,i) = str2double(data(3)) ;
        Z(t,i) = str2double(data(4)) ;
    end

end
fclose(fid) ;

% COMPUTE Q AND Rg
Q = zeros(1,T) ; Rg = zeros(1,T) ;
for t=1:T

    Q(t) = sqrt((X(t,end)-X(t,1))^2+(Y(t,end)-Y(t,1))^2+(Z(t,end)-Z(t,1))^2) ;

    Rcm = [mean(X(t,:)) mean(Y(t,:)) mean(Z(t,:))] ;
    d2 = (X(t,:)-Rcm(1)).^2+(Y(t,:)-Rcm(2)).^2+(Z(t,:)-Rcm(3)).^2 ;
    Rg(t) = sqrt(mean(d2)) ;

end

% PLOT
figure(1) ;
plot(Q,'bo-') ;
xlabel('index') ;
ylabel('end-to-end distance Q') ;

figure(2) ;
plot(Rg, 'ro-') ;
xlabel('index') ;
ylabel('radius of gyration R_g') ;

% MEAN SQUARES
mQ2_sim = mean(Q.^2)
mQ2_th = N*b*b

mRg2_sim = mean(Rg.^2)
mRg2_th = N*b*b/6

%Up until here it's basically the teachers code. The histogram is pretty easy to comptue
fig3 = figure;

histogram(Q, bins, ...
    'Normalization','pdf', ...
    'FaceColor',[0.7 0.7 0.9], ...
    'EdgeColor','k');
hold on;

title('Histogram of end-to-end distance Q');
xlabel('Q');
ylabel('Probability density');

% Theoretical distribution
x = linspace(0, max(Q)*1.1, 400);

P = 4*pi .* x.^2 .* ...
    (3/(2*pi*N*b^2))^(3/2) .* ...
    exp(-(3*x.^2)/(2*N*b^2));

plot(x, P, 'r', 'LineWidth',2);

legend('Simulation','Theory');
grid on;

exportgraphics(fig3, 'Histogram_Q.pdf', 'ContentType','vector');
