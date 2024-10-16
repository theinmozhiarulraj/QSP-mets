figure();

tiledlayout(5,4);
nexttile
simbio_plot(simData,'APC' ,'CompartmentName','V_T' ,'LegendEntry','$APC_{T}$'  );
% xlim([0 200]);

nexttile
simbio_plot(simData,'APC' ,'CompartmentName','V_T_Ln1','LegendEntry','$APC_{T-Ln}$' );
% xlim([0 200]);

nexttile
simbio_plot(simData,'mAPC','CompartmentName','V_T' ,'LegendEntry','$mAPC_{T}$' );
% xlim([0 200]);

nexttile
simbio_plot(simData,'mAPC','CompartmentName','V_T_Ln1','LegendEntry','$mAPC_{T-Ln}$');
% xlim([0 200]);

nexttile
simbio_plot(simData,'T0','CompartmentName','V_T' ,'LegendEntry','$T0_{T}$' );
% xlim([0 200]);

nexttile
simbio_plot(simData,'T0','CompartmentName','V_T_Ln1','LegendEntry','$T0_{T-Ln}$');
% xlim([0 200]);

nexttile
simbio_plot(simData,'T1','CompartmentName','V_T' ,'LegendEntry','$T1_{T}$' );
%xlim([0 200]);

nexttile
simbio_plot(simData,'T1','CompartmentName','V_T_Ln1','LegendEntry','$T1_{T-Ln}$');
%xlim([0 200]);

nexttile
simbio_plot(simData,'Th','CompartmentName','V_T' ,'LegendEntry','$Th_{T}$' );
%xlim([0 200]);

nexttile
simbio_plot(simData,'Th','CompartmentName','V_T_Ln1','LegendEntry','$Th_{T-Ln}$');
%xlim([0 200]);

nexttile
simbio_plot(simData,'Mac_M1','CompartmentName','V_T' ,'LegendEntry','$Mac-M1_{T}$' );
%xlim([0 200]);

nexttile
simbio_plot(simData,'Mac_M1','CompartmentName','V_T_Ln1','LegendEntry','$Mac-M1_{T-Ln}$');
%xlim([0 200]);
