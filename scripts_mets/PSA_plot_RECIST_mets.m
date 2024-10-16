% SimBiology Plot tumor size for PSA cases and RECIST criteria
%
% Inputs: simDataPSA     -- Object containing the simbiology model outputs
%                            for all batch simulations
%         simDataPSApost -- Object containing the post processed simbiology
%                           model outputs for all batch simulations
%         params_out     -- object containing model outputs with RCEIST and
%                           Response status


function PSA_plot_RECIST_mets(simDataPSA,simDataPSApost,params_out)

% n_PSA = length(params_out.RECIST(params_out.iPlausSim));
n_PSA = length(params_out.iPatientPlaus);
index = params_out.iPatientPlaus;

npatients_plotting=100;
perm_index=randperm(npatients_plotting);
perm_patients=index(perm_index);
figure; hold on; box on;

j = [find(strcmp(simDataPSApost(perm_patients(1)).simData.DataNames,'D_T_all_perc') )];

for i =1:npatients_plotting
    plot(simDataPSA(perm_patients(i)).simData.Time/30, ...
         simDataPSApost(perm_patients(i)).simData.Data(:,j) ,'DisplayName',num2str(perm_patients(i)), 'LineWidth', 2)
end
time = simDataPSA(index(1)).simData.Time/30;
plot([0 max(time)*1.25], [+20 +20] , '--k' )
plot([0 max(time)*1.25], [-30 -30] , '--k' )
hx = text(max(time)*1.10, 60,'PD');
hx = text(max(time)*1.10,  0,'SD');
hx = text(max(time)*1.05,-60,'PR/CR');
xlabel('Time (month)','Fontsize',16);
ylabel('% change in sum of tumor diameters','Fontsize',16);
ylim([-100 100])
xlim([0  max(time)*1.25])
set(gca,'Fontsize',14)