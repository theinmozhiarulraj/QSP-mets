% Function to generate object with model parameters to include in parameter
% sensitivity analysis
%
% Output: params -- object containing parameters
%                   -> for each parameter:
%                       - Adds the name of the parameter to the list
%                       - defines upper and lower bounds for uniform and
%                       loguniform
%                       - defines Median and Sigma for normal and lognormal
%                       - specifies the Sampling technique choose from:
%                           - uniform
%                           - loguniform
%                           - normal
%                           - lognormal
%         Examples:
%             % k1
%             params.names = [params.names; 'k1'];
%             params.k1.UpperBound = 1;
%             params.k1.LowerBound = 0;
%             params.k1.Sampling   = 'uniform';
%             params.k1.ScreenName = 'k1 binding rate';
%
%             % k2
%             params.names = [params.names; 'k2'];
%             params.k2.UpperBound = 1;
%             params.k2.LowerBound = 0;
%             params.k2.Sampling   = 'loguniform';
%             params.k2.ScreenName = 'k2 binding rate';
%
%             % k3
%             params.names = [params.names; 'k3'];
%             params.k3.Median     = 1;
%             params.k3.Sigma      = 1;
%             params.k3.Sampling   = 'normal';
%             params.k3.ScreenName = 'k3 binding rate';
%
%             % k4
%             params.names = [params.names; 'k4'];
%             params.k4.Median     = 1;
%             params.k4.Sigma      = 1;
%             params.k4.Sampling   = 'lognormal';
%             params.k4.ScreenName = 'k4 binding rate';

function params = PSA_param_in_TNBC_w60(TumorArray,cancer_clones,n_T_specs)

ntumors=length(TumorArray);
ncancer_clones=length(cancer_clones);

params.names = {};

% add fold changes here
% params.names = [params.names; 'foldchange_k_C1_growth_Ln1'];
% params.foldchange_k_C1_growth_Ln1.UpperBound = 1;
% params.foldchange_k_C1_growth_Ln1.LowerBound = 1;
% params.foldchange_k_C1_growth_Ln1.Sampling   = 'uniform';
% params.foldchange_k_C1_growth_Ln1.ScreenName = 'Fold change in growth rate';

% Initial Tumour Diameter of atleast 1 met
params.names = [params.names; 'initial_met_diameter'];
params.initial_met_diameter.Median = log(1.65);
params.initial_met_diameter.Sigma = 0.3;
params.initial_met_diameter.Sampling   = 'lognormal';
params.initial_met_diameter.ScreenName = 'initial met diameter';

% thein: sample fractions for t cell specificities
T_specs_fractions=sample_fractions(n_T_specs); 

for i=1:ncancer_clones
    for j=1:n_T_specs
        params.names = [params.names; ['agconc_' cancer_clones{i} '_' num2str(j)]];
        params.(['agconc_' cancer_clones{i} '_' num2str(j)]).Median = log(5.4e-13/n_T_specs);
        params.(['agconc_' cancer_clones{i} '_' num2str(j)]).Sigma = 1;
        params.(['agconc_' cancer_clones{i} '_' num2str(j)]).Sampling   = 'lognormal';
        params.(['agconc_' cancer_clones{i} '_' num2str(j)]).ScreenName = 'concentration of neoantigens';
    end
end


% initial number of cancer cells of each clone of each tumor
% 
for j=1:ncancer_clones
    params.names = [params.names; ['ncells_' cancer_clones{j}]];
    params.(['ncells_' cancer_clones{j}]).Median = log(4.7e6/ncancer_clones); %4.7e6*clone_fractions(j)
    params.(['ncells_' cancer_clones{j}]).Sigma = 0.3;
    params.(['ncells_' cancer_clones{j}]).Sampling   = 'lognormal';
    params.(['ncells_' cancer_clones{j}]).ScreenName = 'initial cancer cell number of each clone';

end

for i=2:ntumors
    %clone_fractions=sample_fractions(ncancer_clones); 
    for j=1:ncancer_clones
        params.names = [params.names; ['ncells_' cancer_clones{j} TumorArray{i}]];
        params.(['ncells_' cancer_clones{j} TumorArray{i}]).Median = log(4.7e6/ncancer_clones); %4.7e6*clone_fractions(j)
        params.(['ncells_' cancer_clones{j} TumorArray{i}]).Sigma = 0.3;
        params.(['ncells_' cancer_clones{j} TumorArray{i}]).Sampling   = 'lognormal';
        params.(['ncells_' cancer_clones{j} TumorArray{i}]).ScreenName = 'initial cancer cell number of each clone';
    end
end

% for i=2:ntumors
%     params.names = [params.names; ['delay' TumorArray{i}]];
%     params.(['delay' TumorArray{i}]).UpperBound = 30; %30
%     params.(['delay' TumorArray{i}]).LowerBound = 500; %500
%     params.(['delay' TumorArray{i}]).Sampling   = 'uniform';
%     params.(['delay' TumorArray{i}]).ScreenName = 'Time delay of met seeding';
% end

for i=2:3
    params.names = [params.names; ['delay' TumorArray{i}]];
    params.(['delay' TumorArray{i}]).Median = log(2100); %30
    params.(['delay' TumorArray{i}]).Sigma = 1; %500
    params.(['delay' TumorArray{i}]).Sampling   = 'lognormal';
    params.(['delay' TumorArray{i}]).ScreenName = 'Time delay of met seeding';
end

for i=4:4
    params.names = [params.names; ['delay' TumorArray{i}]];
    params.(['delay' TumorArray{i}]).Median = log(2100); %30
    params.(['delay' TumorArray{i}]).Sigma = 1; %500
    params.(['delay' TumorArray{i}]).Sampling   = 'lognormal';
    params.(['delay' TumorArray{i}]).ScreenName = 'Time delay of met seeding';
end


params.names = [params.names; 'surgery_time'];
params.surgery_time.UpperBound = 1;
params.surgery_time.LowerBound = 400;
params.surgery_time.Sampling   = 'uniform';
params.surgery_time.ScreenName = 'Time delay of surgery';

% Cancer cell growth rate 
for nc=1:ncancer_clones
    p_name=['k_' cancer_clones{nc} '_growth'];
    params.names = [params.names; p_name];
    params.(['k_' cancer_clones{nc} '_growth']).Median = log(0.0065); 
    params.(['k_' cancer_clones{nc} '_growth']).Sigma = 0.7; 
    params.(['k_' cancer_clones{nc} '_growth']).Sampling   = 'lognormal';
    params.(['k_' cancer_clones{nc} '_growth']).ScreenName = 'Rate of tumor growth';
end

for i=2:ntumors
    for nc=1:ncancer_clones
        params.names = [params.names; ['k_' cancer_clones{nc} '_growth' TumorArray{i}]];
        params.(['k_' cancer_clones{nc} '_growth' TumorArray{i}]).Median = log(0.0065); 
        params.(['k_' cancer_clones{nc} '_growth' TumorArray{i}]).Sigma = 0.7; 
        params.(['k_' cancer_clones{nc} '_growth' TumorArray{i}]).Sampling   = 'lognormal';
        params.(['k_' cancer_clones{nc} '_growth' TumorArray{i}]).ScreenName = 'Rate of tumor growth';
    end
end

% % C1 basal death rate (estimated)
% for nc=1:ncancer_clones
%     params.names = [params.names; ['k_' cancer_clones{nc} '_death']];
%     params.(['k_' cancer_clones{nc} '_death']).UpperBound = 0.001; %0.001
%     params.(['k_' cancer_clones{nc} '_death']).LowerBound = 0.00001; % 0.00001
%     params.(['k_' cancer_clones{nc} '_death']).Sampling   = 'loguniform';
%     params.(['k_' cancer_clones{nc} '_death']).ScreenName = 'Rate of apoptotic cancer cell death';
%     
%     for i=2:ntumors
%         params.names = [params.names; ['k_' cancer_clones{nc} '_death' TumorArray{i}]];
%         params.(['k_' cancer_clones{nc} '_death' TumorArray{i}]).UpperBound = 0.001; %0.001
%         params.(['k_' cancer_clones{nc} '_death' TumorArray{i}]).LowerBound = 0.00001; %0.00001
%         params.(['k_' cancer_clones{nc} '_death' TumorArray{i}]).Sampling   = 'loguniform';
%         params.(['k_' cancer_clones{nc} '_death' TumorArray{i}]).ScreenName = 'Rate of apoptotic cancer cell death';
%     end
% end
% 
% % T cell exhaustion by cancer cell (estimated)
% params.names = [params.names; 'k_T1'];
% params.k_T1.UpperBound = 1; 
% params.k_T1.LowerBound = 0.01; 
% params.k_T1.Sampling   = 'loguniform';
% params.k_T1.ScreenName = 'Rate of T cell exhaustion by cancer cells';
% 
% % T cell killing of cancer cell (26872694)

% C1 basal death rate
% params.names = [params.names; 'k_C1_death'];
% params.k_C1_death.UpperBound = 0.001;
% params.k_C1_death.LowerBound = 0.00001;
% params.k_C1_death.Sampling   = 'loguniform';
% params.k_C1_death.ScreenName = 'Rate of apoptotic cancer cell death';

% T cell exhaustion by cancer cell
% params.names = [params.names; 'k_T1'];
% params.k_T1.UpperBound = 1;
% params.k_T1.LowerBound = 0.01;
% params.k_T1.Sampling   = 'loguniform';
% params.k_T1.ScreenName = 'Rate of T cell exhaustion by cancer cells';

% Treg inhibition of Teff
% params.names = [params.names; 'k_Treg'];
% params.k_Treg.UpperBound = 1;
% params.k_Treg.LowerBound = 0.01;
% params.k_Treg.Sampling   = 'loguniform';
% params.k_Treg.ScreenName = 'Rate of Teff inhibition by Treg';

% T cell killing of cancer cell
params.names = [params.names; 'k_C_T1'];
params.k_C_T1.Median = log(0.95);
params.k_C_T1.Sigma = 1; 
params.k_C_T1.Sampling   = 'lognormal';
params.k_C_T1.ScreenName = 'Rate of cancer cell death by Teffs';

% % Kd of antigen P1
% params.names = [params.names; 'k_P1_d1'];
% params.k_P1_d1.Median = log(27e-9);
% params.k_P1_d1.Sigma = 1;
% params.k_P1_d1.Sampling   = 'lognormal';
% params.k_P1_d1.ScreenName = 'Kd of tumor neoantigens';

for t=1:n_T_specs
    params.names = [params.names; ['k_P' num2str(t) '_d1']];
    params.(['k_P' num2str(t) '_d1']).Median = log(27e-9);
    params.(['k_P' num2str(t) '_d1']).Sigma = 1;
    params.(['k_P' num2str(t) '_d1']).Sampling   = 'lognormal';
    params.(['k_P' num2str(t) '_d1']).ScreenName = 'Kd of tumor neoantigens';
end

for t=1:n_T_specs
    params.names = [params.names; ['n_T' num2str(t) '_clones']];
    params.(['n_T' num2str(t) '_clones']).Median = log(63*3); 
    params.(['n_T' num2str(t) '_clones']).Sigma = 0.7;
    params.(['n_T' num2str(t) '_clones']).Sampling   = 'lognormal';
    params.(['n_T' num2str(t) '_clones']).ScreenName = 'Tumor-specific T cell clone';
end

params.names = [params.names; 'n_T0_clones'];
params.n_T0_clones.Median = log(63);
params.n_T0_clones.Sigma = 0.7;
params.n_T0_clones.Sampling   = 'lognormal';
params.n_T0_clones.ScreenName = 'Self-antigen-specific T cell clone';

% Number of T Cell Divisions during Activation by IL-2
params.names = [params.names; 'N_IL2_CD8'];
params.N_IL2_CD8.LowerBound = 10;
params.N_IL2_CD8.UpperBound = 12;
params.N_IL2_CD8.Sampling   = 'uniform';
params.N_IL2_CD8.ScreenName = '# of CD8 T cell division by IL-2';

params.names = [params.names; 'N_IL2_CD4'];
params.N_IL2_CD4.LowerBound = 7;
params.N_IL2_CD4.UpperBound = 10;
params.N_IL2_CD4.Sampling   = 'uniform';
params.N_IL2_CD4.ScreenName = '# of CD4 T cell division by IL-2';

% Initial tumor diameter
params.names = [params.names; 'initial_tumour_diameter'];
params.initial_tumour_diameter.Median = log(2.5);
params.initial_tumour_diameter.Sigma = 0.4;
params.initial_tumour_diameter.Sampling   = 'lognormal';
params.initial_tumour_diameter.ScreenName = 'Initial tumor diameter';

% Checkpoint Module
% Average PD1 Expression on T cells
params.names = [params.names; 'T_PD1_total'];
params.T_PD1_total.UpperBound = 3.1e3*20;
params.T_PD1_total.LowerBound = 1.2e4;
params.T_PD1_total.Sampling   = 'loguniform';
params.T_PD1_total.ScreenName = 'Average number of PD1 per T cell';

% Average Baseline PDL1 Expression on Tumor/Immune cells in tumor
params.names = [params.names; 'C1_PDL1_base'];
params.C1_PDL1_base.UpperBound = 8e4*20/6;
params.C1_PDL1_base.LowerBound = 5e4;
params.C1_PDL1_base.Sampling   = 'loguniform';
params.C1_PDL1_base.ScreenName = 'Baseline number of PDL1 per cell in the tumor';
% Average Baseline PDL1 Expression on mAPCs in TDLN
params.names = [params.names; 'APC_PDL1_base'];
params.APC_PDL1_base.UpperBound = 8e4*20/6;
params.APC_PDL1_base.LowerBound = 5e4;
params.APC_PDL1_base.Sampling   = 'loguniform';
params.APC_PDL1_base.ScreenName = 'Baseline number of PDL1 per mAPC';

% PDL1 Expression on T Cells
% params.names = [params.names; 'T8_PDL1'];
% params.T8_PDL1.UpperBound = 9.3e3*20;
% params.T8_PDL1.LowerBound = 5e4;
% params.T8_PDL1.Sampling   = 'loguniform';
% params.T8_PDL1.ScreenName = 'Baseline number of PDL1 per T cell';

% Effective Density of Bound PD-1 on Teff and Phagocytosis Inhibition
params.names = [params.names; 'PD1_50'];
params.PD1_50.Median = log(6);
params.PD1_50.Sigma = 1;
params.PD1_50.Sampling   = 'lognormal';
params.PD1_50.ScreenName = 'Half-maximal PD1-PDL1 for Teff inhibition';

% Ratio between PDL1 and PDL2 expression in tumor
params.names = [params.names; 'r_PDL2C1'];
params.r_PDL2C1.UpperBound = 0.1;
params.r_PDL2C1.LowerBound = 0.01;
params.r_PDL2C1.Sampling   = 'loguniform';
params.r_PDL2C1.ScreenName = 'PDL2/PDL1 ratio on cancer cell';

% Ratio between PDL1 and PDL2 expression on APC
params.names = [params.names; 'r_PDL2APC'];
params.r_PDL2APC.UpperBound = 0.1;
params.r_PDL2APC.LowerBound = 0.01;
params.r_PDL2APC.Sampling   = 'loguniform';
params.r_PDL2APC.ScreenName = 'PDL2/PDL1 ratio on APCs';

% Maximal rate of Th-to-Treg transdifferentiation
params.names = [params.names; 'k_Th_Treg'];
params.k_Th_Treg.Median = log(0.022);
params.k_Th_Treg.Sigma = 1;
params.k_Th_Treg.Sampling   = 'lognormal';
params.k_Th_Treg.ScreenName = 'Rate of Th to Treg differentiation';
% Rate of tumor vasculature growth fitted to in vivo MX-1 growth
params.names = [params.names; 'k_K_g'];
params.k_K_g.UpperBound = 2.12;
params.k_K_g.LowerBound = 6.12;
params.k_K_g.Sampling   = 'uniform';
params.k_K_g.ScreenName = 'Rate of tumor vasculature growth';

%% Macrophages/MDSC Module
% MDSC Recruitment
params.names = [params.names; 'k_MDSC_mig'];
params.k_MDSC_mig.Median = log(1.1e4);
params.k_MDSC_mig.Sigma = 0.6;
params.k_MDSC_mig.Sampling   = 'lognormal';
params.k_MDSC_mig.ScreenName = 'Recruitment rate of MDSC ';

for i=2:3
    params.names = [params.names; ['k_MDSC_mig' TumorArray{i}]];
    params.(['k_MDSC_mig' TumorArray{i}]).Median = log(1.1e4);
    params.(['k_MDSC_mig' TumorArray{i}]).Sigma = 0.6;
    params.(['k_MDSC_mig' TumorArray{i}]).Sampling   = 'lognormal';
    params.(['k_MDSC_mig' TumorArray{i}]).ScreenName = 'Recruitment rate of MDSC ';
end

for i=4:4
    params.names = [params.names; ['k_MDSC_mig' TumorArray{i}]];
    params.(['k_MDSC_mig' TumorArray{i}]).Median = log(6.5e4);
    params.(['k_MDSC_mig' TumorArray{i}]).Sigma = 0.6;
    params.(['k_MDSC_mig' TumorArray{i}]).Sampling   = 'lognormal';
    params.(['k_MDSC_mig' TumorArray{i}]).ScreenName = 'Recruitment rate of MDSC ';
end

% Macrophage Recruitment
params.names = [params.names; 'k_Mac_mig'];
params.k_Mac_mig.Median = log(1.7e5);
params.k_Mac_mig.Sigma = 0.6;
params.k_Mac_mig.Sampling   = 'lognormal';
params.k_Mac_mig.ScreenName = 'Recruitment rate of macrophage ';

for i=2:3
    params.names = [params.names; ['k_Mac_mig' TumorArray{i}]];
    params.(['k_Mac_mig' TumorArray{i}]).Median = log(3.5e5); 
    params.(['k_Mac_mig' TumorArray{i}]).Sigma = 0.6;
    params.(['k_Mac_mig' TumorArray{i}]).Sampling   = 'lognormal';
    params.(['k_Mac_mig' TumorArray{i}]).ScreenName = 'Recruitment rate of macrophage ';
end

for i=4:4
    params.names = [params.names; ['k_Mac_mig' TumorArray{i}]];
    params.(['k_Mac_mig' TumorArray{i}]).Median = log(3.0e5); % 1.7e5
    params.(['k_Mac_mig' TumorArray{i}]).Sigma = 0.6;
    params.(['k_Mac_mig' TumorArray{i}]).Sampling   = 'lognormal';
    params.(['k_Mac_mig' TumorArray{i}]).ScreenName = 'Recruitment rate of macrophage ';
end

% CCL2 Secretion
params.names = [params.names; 'k_CCL2_sec'];
params.k_CCL2_sec.Median = log(1.7e-12);
params.k_CCL2_sec.Sigma = 0.6;
params.k_CCL2_sec.Sampling   = 'lognormal';
params.k_CCL2_sec.ScreenName = 'Secretion rate of CCL2 ';

% CD47 Expression
params.names = [params.names; 'C_CD47'];
params.C_CD47.UpperBound = 700;
params.C_CD47.LowerBound = 100;
params.C_CD47.Sampling = 'uniform';
params.C_CD47.ScreenName = 'Mean CD47 density on cancer cell ';
% SIRPa Expression
params.names = [params.names; 'M_SIRPa'];
params.M_SIRPa.UpperBound = 180;
params.M_SIRPa.LowerBound = 20;
params.M_SIRPa.Sampling   = 'uniform';
params.M_SIRPa.ScreenName = 'Mean SIRPa density on macrophage ';
% TAM PD1 Expression
params.names = [params.names; 'M_PD1_total'];
params.M_PD1_total.UpperBound = 3.1e3*20;
params.M_PD1_total.LowerBound = 1.5e3;
params.M_PD1_total.Sampling   = 'loguniform';
params.M_PD1_total.ScreenName = 'Mean PD1 density on macrophage ';
% Maximal rate macrophage-mediated phagocytosis
params.names = [params.names; 'k_M1_phago'];
params.k_M1_phago.Median = log(.33);
params.k_M1_phago.Sigma = 1;
params.k_M1_phago.Sampling   = 'lognormal';
params.k_M1_phago.ScreenName = 'Maximal rate of macrophage-mediated phagocytosis ';

params.names = [params.names; 'SIRPa_50'];
params.SIRPa_50.Median = log(37);
params.SIRPa_50.Sigma = .5;
params.SIRPa_50.Sampling   = 'lognormal';
params.SIRPa_50.ScreenName = 'Half-maximal SIRPa-CD47 for phagocytosis inhibition';

params.names = [params.names; 'K_Mac_C'];
params.K_Mac_C.Median = log(2);
params.K_Mac_C.Sigma = 1;
params.K_Mac_C.Sampling   = 'lognormal';
params.K_Mac_C.ScreenName = 'Dependence of phagocytosis on M1/C ratio';


%% Nab-paclitaxel
% nabPaclitaxel PK Parameters
% params.names = [params.names; 'Vmcl'];
% params.Vmcl.UpperBound = 9836;
% params.Vmcl.LowerBound = 6500;
% params.Vmcl.Sampling   = 'uniform';
% params.Vmcl.ScreenName = 'Max clearance rate from V1 compartment ';
% 
% params.names = [params.names; 'Kcl'];
% params.Kcl.UpperBound = 58.9;
% params.Kcl.LowerBound = 24.9;
% params.Kcl.Sampling   = 'uniform';
% params.Kcl.ScreenName = 'Half-max conc. of Nab-P for V1 clearance ';
% 
% params.names = [params.names; 'Vmt'];
% params.Vmt.UpperBound = 540445;
% params.Vmt.LowerBound = 190694;
% params.Vmt.Sampling   = 'uniform';
% params.Vmt.ScreenName = 'Max clearance rate from V1 to V2 distribution ';
% 
% params.names = [params.names; 'Kt'];
% params.Kt.UpperBound = 7910;
% params.Kt.LowerBound = 2210;
% params.Kt.Sampling   = 'uniform';
% params.Kt.ScreenName = 'Half-max conc. of Nab-P for V1 to V2 distribution ';
% 
% params.names = [params.names; 'BSA'];
% params.BSA.UpperBound = 2.4;
% params.BSA.LowerBound = 1.3;
% params.BSA.Sampling   = 'uniform';
% params.BSA.ScreenName = 'Body surface area ';
% 
% params.names = [params.names; 'V_1'];
% params.V_1.UpperBound = 17.85;
% params.V_1.LowerBound = 13.71;
% params.V_1.Sampling   = 'uniform';
% params.V_1.ScreenName = 'Peripheral compartment (V1) for Nab-P PK';
% 
% params.names = [params.names; 'V_2'];
% params.V_2.UpperBound = 1935;
% params.V_2.LowerBound = 1396;
% params.V_2.Sampling   = 'uniform';
% params.V_2.ScreenName = 'Peripheral compartment (V2) for Nab-P PK ';
% 
% params.names = [params.names; 'V_3'];
% params.V_3.UpperBound = 99.1;
% params.V_3.LowerBound = 59.8;
% params.V_3.Sampling   = 'uniform';
% params.V_3.ScreenName = 'Peripheral compartment (V3) for Nab-P PK ';
% 
% % PD parameters
% params.names = [params.names; 'r_nabp'];
% params.r_nabp.UpperBound = 2;
% params.r_nabp.LowerBound = 1;
% params.r_nabp.Sampling   = 'uniform';
% params.r_nabp.ScreenName = 'Tumor to plasma conc. ratio of Nab-P ';
% 
% params.names = [params.names; 'IC50_nabp'];
% params.IC50_nabp.Median = log(92);
% params.IC50_nabp.Sigma = 1;
% params.IC50_nabp.Sampling   = 'lognormal';
% params.IC50_nabp.ScreenName = 'Half-max conc. of Nab-P for cancer killing ';
% 
% params.names = [params.names; 'k_C_resist'];
% params.k_C_resist.Median = log(1e-4);
% params.k_C_resist.Sigma = 1;
% params.k_C_resist.Sampling   = 'lognormal';
% params.k_C_resist.ScreenName = 'Rate of chemo-resistance development ';
% 
% params.names = [params.names; 'k_vas_nabp'];
% params.k_vas_nabp.Median = log(2.8e-4);
% params.k_vas_nabp.Sigma = 1;
% params.k_vas_nabp.Sampling   = 'lognormal';
% params.k_vas_nabp.ScreenName = 'Rate of angiogenic factor induction by Nab-P ';
% 
% params.names = [params.names; 'Kc_nabp'];
% params.Kc_nabp.Median = log(8e7);
% params.Kc_nabp.Sigma = 1;
% params.Kc_nabp.Sampling   = 'lognormal';
% params.Kc_nabp.ScreenName = 'Depth of tumor penetrated by Nab-P ';


end
