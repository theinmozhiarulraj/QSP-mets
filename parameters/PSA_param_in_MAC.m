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

function params = PSA_param_in_MAC

params.names = {};

params.names = [params.names; 'k_Mac_mig'];
params.k_Mac_mig.Median = log(1.7e5);
params.k_Mac_mig.Sigma = 1;
params.k_Mac_mig.Sampling   = 'lognormal';
params.k_Mac_mig.ScreenName = 'M1 macrophage migration';

params.names = [params.names; 'k_Mac_death'];
params.k_Mac_death.Median = log(0.02);
params.k_Mac_death.Sigma = 1;
params.k_Mac_death.Sampling   = 'lognormal';
params.k_Mac_death.ScreenName = 'Macrophage death';

params.names = [params.names; 'k_TGFb_Msec'];
params.k_TGFb_Msec.Median = log(2e-11);
params.k_TGFb_Msec.Sigma = 1;
params.k_TGFb_Msec.Sampling   = 'lognormal';
params.k_TGFb_Msec.ScreenName = 'TGF\beta secretion by M2';

params.names = [params.names; 'k_vas_Msec'];
params.k_vas_Msec.Median = log(1.7e-5);
params.k_vas_Msec.Sigma = 1;
params.k_vas_Msec.Sampling   = 'lognormal';
params.k_vas_Msec.ScreenName = 'Angiogenic factor secretion by M2';

params.names = [params.names; 'k_IL12_sec'];
params.k_IL12_sec.Median = log(8.5e-12);
params.k_IL12_sec.Sigma = 1;
params.k_IL12_sec.Sampling   = 'lognormal';
params.k_IL12_sec.ScreenName = 'IL-12 secretion by mAPC';

params.names = [params.names; 'k_IL12_Msec'];
params.k_IL12_Msec.Median = log(5e-13);
params.k_IL12_Msec.Sigma = 1;
params.k_IL12_Msec.Sampling   = 'lognormal';
params.k_IL12_Msec.ScreenName = 'IL-12 secretion by M1';

params.names = [params.names; 'k_IL12_deg'];
params.k_IL12_deg.Median = log(0.0231);
params.k_IL12_deg.Sigma = 1;
params.k_IL12_deg.Sampling   = 'lognormal';
params.k_IL12_deg.ScreenName = 'IL-12 degradation';

params.names = [params.names; 'k_IL10_sec'];
params.k_IL10_sec.Median = log(3e-12);
params.k_IL10_sec.Sigma = 1;
params.k_IL10_sec.Sampling   = 'lognormal';
params.k_IL10_sec.ScreenName = 'IL-10 secretion by M2';

params.names = [params.names; 'k_IL10_deg'];
params.k_IL10_deg.Median = log(4);
params.k_IL10_deg.Sigma = 1;
params.k_IL10_deg.Sampling   = 'lognormal';
params.k_IL10_deg.ScreenName = 'IL-10 degradation';

params.names = [params.names; 'k_M2_pol'];
params.k_M2_pol.Median = log(0.25);
params.k_M2_pol.Sigma = 1;
params.k_M2_pol.Sampling   = 'lognormal';
params.k_M2_pol.ScreenName = 'M1-to-M2 polarization';

params.names = [params.names; 'k_M1_pol'];
params.k_M1_pol.Median = log(0.02);
params.k_M1_pol.Sigma = 1;
params.k_M1_pol.Sampling   = 'lognormal';
params.k_M1_pol.ScreenName = 'M2-to-M1 polarization';

params.names = [params.names; 'IL10_50'];
params.IL10_50.Median = log(8);
params.IL10_50.Sigma = 1;
params.IL10_50.Sampling   = 'lognormal';
params.IL10_50.ScreenName = 'EC50 of IL-10 on polarization';

params.names = [params.names; 'IL12_50'];
params.IL12_50.Median = log(0.14);
params.IL12_50.Sigma = 1;
params.IL12_50.Sampling   = 'lognormal';
params.IL12_50.ScreenName = 'EC50 of IL-12';

params.names = [params.names; 'IFNg_50'];
params.IFNg_50.Median = log(2.9);
params.IFNg_50.Sigma = 1;
params.IFNg_50.Sampling   = 'lognormal';
params.IFNg_50.ScreenName = 'EC50 of IFN\gamma';

params.names = [params.names; 'k_M1_phago'];
params.k_M1_phago.Median = log(.3);
params.k_M1_phago.Sigma = 1;
params.k_M1_phago.Sampling   = 'lognormal';
params.k_M1_phago.ScreenName = 'Phagocytosis by M1';

params.names = [params.names; 'TGFb_50'];
params.TGFb_50.Median = log(0.14);
params.TGFb_50.Sigma = 1;
params.TGFb_50.Sampling   = 'lognormal';
params.TGFb_50.ScreenName = 'EC50 of TGF\beta';

params.names = [params.names; 'CCL2_50'];
params.CCL2_50.Median = log(0.23);
params.CCL2_50.Sigma = 1;
params.CCL2_50.Sampling   = 'lognormal';
params.CCL2_50.ScreenName = 'EC50 of CCL2/MCP-1';

params.names = [params.names; 'IL10_50_phago'];
params.IL10_50_phago.Median = log(270);
params.IL10_50_phago.Sigma = 1;
params.IL10_50_phago.Sampling   = 'lognormal';
params.IL10_50_phago.ScreenName = 'EC50 of IL-10 on phagocytosis';

params.names = [params.names; 'K_Mac_C'];
params.K_Mac_C.Median = log(2);
params.K_Mac_C.Sigma = 1;
params.K_Mac_C.Sampling   = 'lognormal';
params.K_Mac_C.ScreenName = 'Dependence of phagocytosis on M1/C ratio';

end
