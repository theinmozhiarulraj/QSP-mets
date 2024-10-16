%% Immune Oncology Model
% Script for setting up and running the immune oncology model in simbiology
clear
%close all
sbioreset

%% Add use-defined units into both Simbiology and Symbolic libraries
% Add 'cell' unit to SimBiology and Symbolic Toolboxes
if (isempty(sbioshowunits('cell')))
    cell_unit = sbiounit('cell','molecule');
    sbioaddtolibrary(cell_unit);
end
% Symbolic Unit
u = symunit;
try u.cell;
catch
    newUnit('cell',u.molecule);
end

% Add 'mU' unit to SimBiology and Symbolic Toolboxes
mU_unit = sbiounit('mU','mole/liter');
sbioaddtolibrary(mU_unit);
% Symbolic Unit
u = symunit;
try u.mU;
catch
    newUnit('mU',u.molarity);
end

%  thein: only primary tumor
%  TumorArray = {''}; 
%  LNArray = {''};

%  thein: primary +1 lung met.
%  TumorArray = {'','_Ln1'}; 
%  LNArray = {'','_Ln1'}; 

% thein: primary + 2 lung mets
%  TumorArray = {'','_Ln1','_Ln2'}; 
%  LNArray = {'','_Ln1','_Ln1'}; 

%  thein: primary + 2 lung met. + 1 other
TumorArray = {'','_Ln1','_Ln2','_other'}; 
LNArray = {'','_Ln1','_Ln1','_Ln2'}; 

% thein: number of cancer clones
ncancer_clones=5;
cancer_clones={};
for nc=1:ncancer_clones
    cancer_clones=[cancer_clones,['C' num2str(nc)]];
end

% thein: number of T cell specificities/neoantigens
n_T_specs=8; 

% thein: sample Ag proportions
%rng(1)
% ag_fraction_clones = sample_ag_proportions(ncancer_clones,n_T_specs);

%% thein: use this to automate the met naming

% number of lung mets
% nlungmets = 4;

% first element of TumorArray and LNArray correspond to primary tumor

% TumorArray = {''};
% for i=1:nlungmets
%     str=['_Ln' num2str(i,'%d')];
%     TumorArray = [TumorArray; str];
% end
% TumorArray;
% 
% LNArray = repelem([{''}, {'_Ln1'}], [1 nlungmets])';
% 
% if(length(TumorArray)~=length(LNArray))
%     error('error: length of Tumor and LN arrays should be same');
% end

%% Setup Parameters
% Setup Parameters
params_in     = parameters_TNBC_mets_w60_vct(TumorArray,cancer_clones,n_T_specs);
params_out    = load_parameters(params_in);

%% Create the SimBiology Model
% Model Settings
model_name = 'Immune Oncology Model';
start_time = 0.0; % [days]
time_step = 1; % [days] 0.01 days ~ 15 mins
%end_time = 400; % [days]
end_time = 400; % [days] 
absolute_tolerance = 1e-9;%1e-9
relative_tolerance = 1e-6;%1e-6
%solver = 'ode15s';
solver = 'sundials';
% Model Object
time = start_time:time_step:end_time;

%model = simbio_init(model_name,time,solver,absolute_tolerance,relative_tolerance,params_out,ncancer_clones,n_T_specs,ag_fraction_clones,TumorArray,LNArray);
model = simbio_init(model_name,time,solver,absolute_tolerance,relative_tolerance,params_out,ncancer_clones,n_T_specs,TumorArray,LNArray);

% Maximal simulation time
config = getconfigset(model);
set(config, 'MaximumWallClock', 120)
%set(config, 'MaximumNumberOfLogs', 1)
set(config.SolverOptions, 'AbsoluteToleranceScaling', false)

%thein: create an array for initial number of cancer cells
ntumors=length(TumorArray);
% ncancer_cells=[];
% for nc=1:ncancer_clones
%     ncells=zeros(1,ntumors);
%     for t=1:ntumors
%        ncells(t)=1e-6;
%        if(t==1)
%           ncells(t)=4.7e6/ncancer_clones; % no primary tumor
%        end
%     end
%     ncancer_cells=[ncancer_cells;ncells];
% end
% disp(ncancer_cells)


%% Add Modules to the Model
% Cancer Modules
model = cancer_module(model,cancer_clones,params_out,TumorArray,LNArray); % 4.7e6

% model = cancer_module(model,'C1',params_out,TumorArray,LNArray,4.7e6); % 4.7e6
% model = cancer_module(model,'C2',params_out,TumorArray,LNArray,{0,0,0,0,0}); % chemotherapy-resistant clone

% T cell Modules
model = Treg_module(model,params_out,n_T_specs,TumorArray,LNArray);
model = Teff_module(model,'1',params_out,cancer_clones,n_T_specs,TumorArray,LNArray);

% APC Module
model = APC_module(model,params_out,TumorArray,LNArray);

% Antigen Modules
%antigenCP = create_antigen(cancer_clones,repelem(5.4e-13,ncancer_clones),'antigenID',0);
antigenCP = create_antigen(model,cancer_clones,'antigenID',0);
model = antigen_module(model,'0',params_out,antigenCP,TumorArray,LNArray);

% antigen   = create_antigen(cancer_clones,[5.4e-13 5.4e-13],'antigenID',1);
% model = antigen_module(model,'1',params_out,antigen,TumorArray,LNArray);

for nt=1:n_T_specs
   %ag_conc = ag_fraction_clones(:,nt)*5.4e-13;
   %antigen = create_antigen(cancer_clones,ag_conc,'antigenID',nt);
   antigen = create_antigen(model,cancer_clones,'antigenID',nt);

   %antigen = create_antigen(cancer_clones,[5.4e-13 5.4e-13],'antigenID',nt);
   model = antigen_module(model,num2str(nt),params_out,antigen,TumorArray,LNArray);
end

% Checkpoint Modules
model = checkpoint_module(model,params_out,'T','C1',TumorArray,LNArray);
model = checkpoint_module(model,params_out,'T','APC',TumorArray,LNArray);
% ADCC Module (use in ipilimumab therapy; in development)
% model = Treg_ADCC_module(model,params_out,TumorArray,LNArray);

% QSPIO-TNBC Modules
model = Th_module(model,params_out,TumorArray,LNArray);

model = MDSC_module(model,params_out,TumorArray,LNArray,cancer_clones,'inostat',0,'drugName','entinostat');

% thein: this module gives an error if there are no cancer cells (C_total=0)
% model = nabpaclitaxel_module(model,params_out,TumorArray,LNArray);
model = macrophage_module(model,params_out,cancer_clones,TumorArray,LNArray,'aCD47',0); % PK module in development for aCD47

%% Setup Dosing
dose_schedule = [];
% dose_schedule = schedule_dosing({'atezolizumab'});
% dose_schedule = schedule_dosing({'nabPaclitaxel'});

% dbstop if warning
simData = sbiosimulate(model,[],[],dose_schedule);

%% Initialize and Run the Model (should run with realistic baseline parameters)
% (should be commented out when conducting in silico virtual clinical trial)
% tic
% [model,success,simDataInit] = initial_conditions(model,TumorArray);
% [model,success] = initial_conditions(model,TumorArray);
% toc
% % 
% % % Generate a list of parameters and species for debug
% % modelComp = listModelComp(model);
% % 
% %model=1;
% % Run Simulation
% if (success)
%     tic
%     simData = sbiosimulate(model,[],[],dose_schedule);
%     toca
% else
%     simData = simDataInit;
%     disp('Tumour did not reach specified initial tumour diameter with current parameters');
% end

%% Plots
% Plot diagnostics
% if (success)
%    diagnostic_plot(simData,model);
%    diagnostic_plot_H(simData);
%    diagnostic_plot_KPR(simData,model);
% end


%  plot(simData.time,simData.Data(:,407),'-k') 
%   hold on;
%  plot(simData.time,simData.Data(:,408),'-g') 
% 
% hold on;
%  simbio_plot(simData,'V_T')
%  hold on;
%  simbio_plot(simData,'V_T_Ln1')
% 
 %figure()
 %simbio_plot(simData,'T_total')
 %hold on;
%  simbio_plot(simData,'C2')
%  hold on;
%  simbio_plot(simData,'C3')

% simbio_plot(simData,'R_Tcell')
% hold on;
% figure()
% simbio_plot(simData,'T1' ,'CompartmentName','V_T' ,'LegendEntry','$APC_{T}$'  );
% hold on;
% simbio_plot(simData,'T2' ,'CompartmentName','V_T' ,'LegendEntry','$APC_{T}$'  );
% hold on;
% simbio_plot(simData,'T3' ,'CompartmentName','V_T' ,'LegendEntry','$APC_{T}$'  );
% hold on;
% simbio_plot(simData,'T4' ,'CompartmentName','V_T' ,'LegendEntry','$APC_{T}$'  );
% hold on;
% simbio_plot(simData,'T5' ,'CompartmentName','V_T' ,'LegendEntry','$APC_{T}$'  );
% hold on;
% 


