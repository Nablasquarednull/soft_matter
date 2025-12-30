
clear all; close all; clc;
% ANALYSIS OF SIMULATION OF 3-D FJC MODEL OF POLYMERS

% PARAMETERS
b = 3.0;
N_values = [10 20 50 100 200 500 700 1000];
T = 1000;

% OUTPUT ARRAYS
mQ2_sim_arr  = zeros(size(N_values));
mRg2_sim_arr = zeros(size(N_values));

for k = 1:length(N_values)

    N = N_values(k);
    fprintf("N = %d...\n", N);

    % INITIALIZATION
    X = zeros(T, N+1);
    Y = zeros(T, N+1);
    Z = zeros(T, N+1);

    % READ INPUT FILE
    filename = sprintf('simulation_FJC_b=%.1f_N=%d_T=%d.xyz', b, N, T);
    fid = fopen(filename, 'r');

    for t = 1:T
        fgetl(fid); % ignore
        fgetl(fid); % ignore

        for i = 1:N+1
            line = fgetl(fid);
            data = strsplit(line);
            X(t,i) = str2double(data{2});
            Y(t,i) = str2double(data{3});
            Z(t,i) = str2double(data{4});
        end
    end
    fclose(fid);

    % COMPUTE Q AND Rg
    Q  = zeros(1,T);
    Rg = zeros(1,T);

    for t = 1:T
        Q(t) = sqrt( (X(t,end)-X(t,1))^2 + ...
                     (Y(t,end)-Y(t,1))^2 + ...
                     (Z(t,end)-Z(t,1))^2 );

        Rcm = [mean(X(t,:)) mean(Y(t,:)) mean(Z(t,:))];
        d2  = (X(t,:)-Rcm(1)).^2 + ...
              (Y(t,:)-Rcm(2)).^2 + ...
              (Z(t,:)-Rcm(3)).^2;

        Rg(t) = sqrt(mean(d2));
    end

    % SAVE MEAN SQUARES
    mQ2_sim_arr(k)  = mean(Q.^2);
    mRg2_sim_arr(k) = mean(Rg.^2);

end
% THEORETICAL SCALING (FJC)
Q2_th  = N_values * b^2;
Rg2_th = N_values * b^2 / 6;

fig1 = figure('Units','inches','Position',[1 1 3.3 3]);

loglog(N_values, mQ2_sim_arr, 'o', ...
       'MarkerSize', 6, 'LineWidth',1.4);
hold on;

loglog(N_values, Q2_th, '--', ...
       'LineWidth',1.6);

xlabel('N');
ylabel('\langle Q^2 \rangle');
legend('Simulation','Theory','Location','northwest');
grid on;

set(gca,'FontName','Times','FontSize',11);
set(gcf,'PaperPositionMode','auto');

exportgraphics(fig1, 'Q2_vs_N_theory.pdf', 'ContentType','vector');

fig2 = figure('Units','inches','Position',[1 1 3.3 3]);

loglog(N_values, mRg2_sim_arr, 'o', ...
       'MarkerSize', 6, 'LineWidth',1.4);
hold on;

loglog(N_values, Rg2_th, '--', ...
       'LineWidth',1.6);

xlabel('N');
ylabel('\langle R_g^2 \rangle');
legend('Simulation','Theory','Location','northwest');
grid on;

set(gca,'FontName','Times','FontSize',11);
set(gcf,'PaperPositionMode','auto');

exportgraphics(fig2, 'Rg2_vs_N_theory.pdf', 'ContentType','vector');

