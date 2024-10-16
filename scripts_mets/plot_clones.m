bar(ag_fraction_clones,'stacked')


for i=1:n_T_specs
    %figure();
    simbio_plot(simData,['T' num2str(i)],'CompartmentName','V_T' ,'LegendEntry','$APC_{T}$'  );
    hold on;
end

for i=1:ncancer_clones
    %figure();
    simbio_plot(simData,['C' num2str(i)],'CompartmentName','V_T' ,'LegendEntry','$APC_{T}$'  );
    hold on;
end
