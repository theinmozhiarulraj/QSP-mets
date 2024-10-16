% load the data from old and new versions
simdata_hanwen = load('matlab_hanwenversion.mat').simData;
simdata_new = load('matlab_newversion.mat').simData;

% number of species
% there are more species in the new version 
[~,nspecies]=size(simdata_hanwen.Data);

nper_col=3;
nper_row=3;
nper_fig=nper_col*nper_row;

start=1;
for j=1:(nspecies/nper_fig)+1

    figure()
    tiledlayout(nper_row,nper_col);

    for i=start:start+nper_fig-1

        start=start+nper_fig; 
    
        compartment=simdata_hanwen.DataInfo{i,1}.Compartment;
        species=simdata_hanwen.DataInfo{i,1}.Name;
    
        nexttile
        simbio_plot(simdata_hanwen,species,'compartmentName',compartment,'legendEntry','Hanwen version');
        hold on;
        simbio_plot(simdata_new,species,'compartmentName',compartment,'legendEntry','New version');
    
        title([species ' ' compartment])
        legend

    end

end


% figure()
% 
% simbio_plot(simdata_hanwen,'P0','compartmentName','V_T','legendEntry','Hanwen version');
% hold on;
% simbio_plot(simdata_new,'P0','compartmentName','V_T','legendEntry','New version');
% xlim([0 10])
% 
% 
% 
% figure()
% 
% simbio_plot(simdata_hanwen,'V_T','legendEntry','Hanwen version');
% hold on;
% simbio_plot(simdata_new,'V_T','legendEntry','New version');
% legend

% load the data 
simdata_prim = load('simData_prim_d100_step1-apc.mat').simData;
simdata_lung = load('simData_lung_d100_step1.mat').simData;
simdata_others = load('simData_others_d100_step1.mat').simData;

delay=100;


figure()
plot(simdata_prim.Time,simdata_prim.Data(:,559),'k')
hold on;
plot(simdata_lung.Time-delay,simdata_lung.Data(:,560),'-g')
hold on;
plot(simdata_others.Time-delay,simdata_others.Data(:,561),':r')
title('Tumor volume')
legend('primary','lung met','other met')

figure()
plot(simdata_prim.Time,simdata_prim.Data(:,425),'k')
hold on;
plot(simdata_lung.Time-delay,simdata_lung.Data(:,436),'-g')
hold on;
plot(simdata_others.Time-delay,simdata_others.Data(:,447),':r')
legend('primary','lung met','other met')

figure()
plot(simdata_prim.Time,simdata_prim.Data(:,428),'k')
hold on;
plot(simdata_lung.Time-delay,simdata_lung.Data(:,439),'-g')
hold on;
plot(simdata_others.Time-delay,simdata_others.Data(:,450),':r')
legend('primary','lung met','other met')

figure()
plot(simdata_prim.Time,simdata_prim.Data(:,427),'k')
hold on;
plot(simdata_lung.Time-delay,simdata_lung.Data(:,438),'-g')
hold on;
plot(simdata_others.Time-delay,simdata_others.Data(:,449),':r')
legend('primary','lung met','other met')

%%%%%%%%% multiple clones

figure()
plot(simdata_prim.Time,simdata_prim.Data(:,890),'k')
hold on;
plot(simdata_lung.Time-delay,simdata_lung.Data(:,891),'-g')
hold on;
plot(simdata_others.Time-delay,simdata_others.Data(:,892),':r')
title('Tumor volume')
legend('primary','lung met','other met')



figure()
plot(simdata_prim.Time,simdata_prim.Data(:,666),'k')
hold on;
plot(simdata_lung.Time-delay,simdata_lung.Data(:,681),'-g')
hold on;
plot(simdata_others.Time-delay,simdata_others.Data(:,696),':r')
legend('primary','lung met','other met')


figure()
plot(simdata_prim.Time,simdata_prim.Data(:,673),'k')
hold on;
plot(simdata_lung.Time-delay,simdata_lung.Data(:,688),'-g')
hold on;
plot(simdata_others.Time-delay,simdata_others.Data(:,703),':r')
legend('primary','lung met','other met')
