%%
clear
close all

immune_oncology_model_TNBC

%%
dose_schedule = [];

n_PSA = 500;

params_in  = PSA_param_in_MAC;
params_out = PSA_param_out(model);
params_in  = PSA_param_obs(params_in);
params_in  = PSA_setup(model,params_in,n_PSA);

warning('off','all')

sbioaccelerate(model, dose_schedule)
tic
[simDataPSA, params_out] = simbio_PSA(model,params_in,params_out,dose_schedule);
toc

% Postprocess Data -> Calculate Clonality, Percentages and ...
simDataPSApost = PSA_post(simDataPSA,params_in,params_out);

% Add pre-treatment observables to the params_in
params_in = PSA_preObs(simDataPSA,simDataPSApost,params_in,params_out);

% Prepare the data for the rest of the analysis
params_out = PSA_prep(simDataPSA,simDataPSApost,params_out);

filename = 'Sensitivity_Analysis_Macrophage.mat';
save(filename, 'params_in', 'params_out', 'simDataPSA', 'simDataPSApost')


%%
n_PSA = length(params_out.iPatientPlaus);
index = params_out.iPatientPlaus;

for i = 1:length(params_in.names)
    namesIn{i} = params_in.(params_in.names{i}).ScreenName;
end
% namesIn{11} = 'M2-to-M1 polarization';
% namesIn{19} = 'Dependence of phagocytosis on M1/C ratio';
% namesIn{18} = 'EC50 of IL-10 on phagocytosis';
% namesIn{17} = 'EC50 of CCL2/MCP-1';


% [rho, pval] = partialcorr(params_in.all(index,:), params_out.post(index,1) , ones(n_PSA,1) , 'rows','complete','type','Spearman');
[rho, pval,~] = PRCC(params_in.all(index,:), params_out.post(index,1), 1, namesIn, 0.05/size(params_in.all(index,:), 2));
[rho_sorted, idx_sorted] = sort(rho);
namesIn_sorted = namesIn(idx_sorted);

figure
hy = barh(rho_sorted, 'FaceColor', 'k');

X = hy.XData;
Y = hy.YData;

temp = ones(size(X,2),1);
temp(Y<0) = -1;

yticks(1:size(X,2))
set(gca,'yticklabel',namesIn_sorted)

% set(gca,'YTickLabelRotation',35)
% set(gca,'XTickLabelRotation',35)
xlim([-.8 .8])
xticks([-.7 -.5 -.3 -.1 .1 .3 .5 .7])
set(gca,'FontSize',14);

%%
for i = 1:size(X,2)
    if pval(i) <= 0.001
        text(Y(i)+.02*temp(i), X(i)-0.3*hy.BarWidth, '***', 'Rotation',90);
    elseif pval(i) <= 0.01
        text(Y(i)+.02*temp(i), X(i)-0.2*hy.BarWidth, '**', 'Rotation',90);
    elseif pval(i) <= 0.05
        text(Y(i)+.02*temp(i), X(i)-0.1*hy.BarWidth, '*', 'Rotation',90);
    end
end

%% phagocytosis module

dose_schedule = [];

n_PSA = 500;

params_in  = PSA_param_in_PHA;
params_out = PSA_param_out(model);
params_in  = PSA_param_obs(params_in);
params_in  = PSA_setup(model,params_in,n_PSA);

warning('off','all')

sbioaccelerate(model, dose_schedule)
tic
[simDataPSA, params_out] = simbio_PSA(model,params_in,params_out,dose_schedule);
toc

% Postprocess Data -> Calculate Clonality, Percentages and ...
simDataPSApost = PSA_post(simDataPSA,params_in,params_out);

% Add pre-treatment observables to the params_in
params_in = PSA_preObs(simDataPSA,simDataPSApost,params_in,params_out);

% Prepare the data for the rest of the analysis
params_out = PSA_prep(simDataPSA,simDataPSApost,params_out);

filename = 'Sensitivity_Analysis_Phagocytosis.mat';
save(filename, 'params_in', 'params_out', 'simDataPSA', 'simDataPSApost')


%%
n_PSA = length(params_out.iPatientPlaus);
index = params_out.iPatientPlaus;

clear H_Mac_C
for i = 1:n_PSA
    [~,temp,~] = selectbyname(simDataPSA(index(i)).simData, 'H_Mac_C');
    H_Mac_C(i,:) = temp;
end

for i = 1:length(params_in.names)
    namesIn{i} = params_in.(params_in.names{i}).ScreenName;
end
% namesIn{1} = 'CD47-SIRP\alpha k_{on}';
% namesIn{2} = 'CD47-SIRP\alpha k_{off}';

% [rho, pval] = partialcorr(params_in.all(index,:), H_Mac_C(index,end) , ones(n_PSA,1) , 'rows','complete','type','Spearman');
[rho, pval,~] = PRCC(params_in.all(index,:), H_Mac_C(index,end), 1, namesIn, 0.05/size(params_in.all(index,:), 2));
[rho_sorted, idx_sorted] = sort(rho);
namesIn_sorted = namesIn(idx_sorted);

figure

hy = barh(rho_sorted, 'FaceColor', 'k');

X = hy.XData;
Y = hy.YData;

temp = ones(size(X,2),1);
temp(Y<0) = -1;

yticks(1:size(X,2))
set(gca,'yticklabel',namesIn_sorted)

% set(gca,'YTickLabelRotation',35)
% set(gca,'XTickLabelRotation',35)
xlim([-.7 .7])
xticks([-.7 -.5 -.3 -.1 .1 .3 .5 .7])
set(gca,'FontSize',14);

%%
for i = 1:size(X,2)
    if pval(i) <= 0.001
        text(Y(i)+.02*temp(i), X(i)-0.3*hy.BarWidth, '***', 'Rotation',90);
    elseif pval(i) <= 0.01
        text(Y(i)+.02*temp(i), X(i)-0.2*hy.BarWidth, '**', 'Rotation',90);
    elseif pval(i) <= 0.05
        text(Y(i)+.02*temp(i), X(i)-0.1*hy.BarWidth, '*', 'Rotation',90);
    end
end
