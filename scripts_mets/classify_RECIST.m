
% calculate percent patients with PD, SD, CR, PR and PR/CR
% only plausible patients are considered

npatients_plausible = length(params_out.iPlaus);

%calc_CR=0;
%calc_PR=0;
calc_PD=0;
calc_CRPR=0;
calc_SD=0;
calc_NP=0; % not patient

for i=1:length(params_out.iPlaus)

    j=params_out.iPlaus(i);

    disp(i)
    if(params_out.RECIST{j}=="CR/PR")
        calc_CRPR=calc_CRPR+1;
    elseif(params_out.RECIST{j}=="SD")
        calc_SD=calc_SD+1;
    elseif(params_out.RECIST{j}=="PD")
        calc_PD=calc_PD+1;
    else
        calc_NP=calc_NP+1;
    end

end

calc_PD=calc_PD/npatients_plausible;
calc_CRPR=calc_CRPR/npatients_plausible;
calc_SD=calc_SD/npatients_plausible;

disp(calc_CRPR)
disp(calc_SD)
disp(calc_PD)



