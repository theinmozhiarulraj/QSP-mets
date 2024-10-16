% Function to generate object with model parameters to include in parameter
% sensitivity analysis
%
% Output: params -- object containing parameters
%                   -> for each parameter:
%                       - Adds the name of the parameter to the list
%                       - defines upper and lower bounds for uniform and
%                       lognormal
%                       - defines Median and Sigma for normal and lognormal
%                       - specifies the Sampling technique choose from:
%                           - uniform
%                           - lognormal
%                           - normal
%                           - lognormal

function params = PSA_param_in_PHA

params.names = {};

% Association rate of CD47-SIRPa
params.names = [params.names; 'kon_CD47_SIRPa'];
params.kon_CD47_SIRPa.Median = log(0.0625);
params.kon_CD47_SIRPa.Sigma = 1;
params.kon_CD47_SIRPa.Sampling = 'lognormal';
params.kon_CD47_SIRPa.ScreenName = 'CD47-SIRP\alpha k_{on}';
% Dissociation rate of CD47-SIRPa
params.names = [params.names; 'koff_CD47_SIRPa'];
params.koff_CD47_SIRPa.Median = log(0.3);
params.koff_CD47_SIRPa.Sigma = 1;
params.koff_CD47_SIRPa.Sampling   = 'lognormal';
params.koff_CD47_SIRPa.ScreenName = 'CD47-SIRP\alpha k_{off}';

% Association rate of CD80-PDL1
params.names = [params.names; 'kon_CD80_PDL1'];
params.kon_CD80_PDL1.Median = log(1.053);
params.kon_CD80_PDL1.Sigma = 1;
params.kon_CD80_PDL1.Sampling = 'lognormal';
params.kon_CD80_PDL1.ScreenName = 'CD80-PDL1 k_{on}';
% Dissociation rate of CD80-PDL1
params.names = [params.names; 'koff_CD80_PDL1'];
params.koff_CD80_PDL1.Median = log(6);
params.koff_CD80_PDL1.Sigma = 1;
params.koff_CD80_PDL1.Sampling   = 'lognormal';
params.koff_CD80_PDL1.ScreenName = 'CD80-PDL1 k_{off}';
% Association rate of PD1-PDL1
params.names = [params.names; 'kon_PD1_PDL1'];
params.kon_PD1_PDL1.Median = log(0.0583);
params.kon_PD1_PDL1.Sigma = 1;
params.kon_PD1_PDL1.Sampling = 'lognormal';
params.kon_PD1_PDL1.ScreenName = 'PD1-PDL1 k_{on}';
% Dissociation rate of PD1-PDL1
params.names = [params.names; 'koff_PD1_PDL1'];
params.koff_PD1_PDL1.Median = log(1.435);
params.koff_PD1_PDL1.Sigma = 1;
params.koff_PD1_PDL1.Sampling   = 'lognormal';
params.koff_PD1_PDL1.ScreenName = 'PD1-PDL1 k_{off}';
% Association rate of PD1-PDL2
params.names = [params.names; 'kon_PD1_PDL2'];
params.kon_PD1_PDL2.Median = log(0.0767);
params.kon_PD1_PDL2.Sigma = 1;
params.kon_PD1_PDL2.Sampling = 'lognormal';
params.kon_PD1_PDL2.ScreenName = 'PD1-PDL2 k_{on}';
% Dissociation rate of PD1-PDL2
params.names = [params.names; 'koff_PD1_PDL2'];
params.koff_PD1_PDL2.Median = log(0.529);
params.koff_PD1_PDL2.Sigma = 1;
params.koff_PD1_PDL2.Sampling   = 'lognormal';
params.koff_PD1_PDL2.ScreenName = 'PD1-PDL2 k_{off}';


% Half-maximal SIRPa binding for phagocytosis inhibition
params.names = [params.names; 'SIRPa_50'];
params.SIRPa_50.Median = log(27);
params.SIRPa_50.Sigma = 1;
params.SIRPa_50.Sampling   = 'lognormal';
params.SIRPa_50.ScreenName = 'EC50 of bound SIRP\alpha';
% Hill coefficient of H_{SIRP\alpha}
params.names = [params.names; 'n_SIRPa'];
params.n_SIRPa.Median = log(1.5);
params.n_SIRPa.Sigma = 1;
params.n_SIRPa.Sampling   = 'lognormal';
params.n_SIRPa.ScreenName = 'Hill coefficient of H_{SIRP\alpha}';
% Half-maximal PD-1 binding for phagocytosis inhibition
params.names = [params.names; 'PD1_50'];
params.PD1_50.Median = log(6);
params.PD1_50.Sigma = 1;
params.PD1_50.Sampling   = 'lognormal';
params.PD1_50.ScreenName = 'EC50 of bound PD-1';
% Hill coefficient of H_{PD1}
params.names = [params.names; 'n_PD1'];
params.n_PD1.Median = log(2);
params.n_PD1.Sigma = 1;
params.n_PD1.Sampling   = 'lognormal';
params.n_PD1.ScreenName = 'Hill coefficient of H_{PD1}';


% CD47 Expression
params.names = [params.names; 'C_CD47'];
params.C_CD47.UpperBound = 700;
params.C_CD47.LowerBound = 100;
params.C_CD47.Sampling = 'uniform';
params.C_CD47.ScreenName = 'Mean CD47 density on C';
% SIRPa Expression
params.names = [params.names; 'M_SIRPa'];
params.M_SIRPa.UpperBound = 180;
params.M_SIRPa.LowerBound = 20;
params.M_SIRPa.Sampling   = 'uniform';
params.M_SIRPa.ScreenName = 'Mean SIRPa density on M';
% TAM PD1 Expression
params.names = [params.names; 'M_PD1_total'];
params.M_PD1_total.UpperBound = 3.1e3*20;
params.M_PD1_total.LowerBound = 1.5e3;
params.M_PD1_total.Sampling   = 'loguniform';
params.M_PD1_total.ScreenName = 'Mean PD1 density on M';
% Average Baseline PDL1 Expression on Tumor/Immune cells in tumor
params.names = [params.names; 'C1_PDL1_base'];
params.C1_PDL1_base.UpperBound = 8e4*20/6; % 5.4e4*20/6 ; 8e4*20/6
params.C1_PDL1_base.LowerBound = 5e4; % 9e3 ; 1e2
params.C1_PDL1_base.Sampling   = 'loguniform';
params.C1_PDL1_base.ScreenName = 'Mean PD-L1 density on C';
% Ratio between PDL1 and PDL2 expression in tumor
params.names = [params.names; 'r_PDL2C1'];
params.r_PDL2C1.UpperBound = 0.1;
params.r_PDL2C1.LowerBound = 0.001;
params.r_PDL2C1.Sampling   = 'loguniform';
params.r_PDL2C1.ScreenName = 'PDL2/PDL1 ratio on C';

end
