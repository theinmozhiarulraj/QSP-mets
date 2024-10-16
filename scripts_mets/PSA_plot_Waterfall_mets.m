% Produce waterfall plot for percent change in tumor size for a parameter
%
% Inputs: simDataPSApost -- Object containing the post processed simbiology
%                           model outputs for all batch simulations
%         params_in      -- Object containing the input parameters
%         param          -- name of the parameter used for color coding


function PSA_plot_Waterfall_mets(simDataPSApost,model,params_in,params_out,param,TumorArray)

test_freq=9*7; %days

n_PSA = length(params_out.iPatientPlaus);
index = params_out.iPatientPlaus;

ntumors=length(TumorArray);
tumor_response=NaN(n_PSA,ntumors);
tumor_pos=zeros(ntumors,1);

for k=1:ntumors
    tumorperc_pos(k,1) = [find(strcmp(simDataPSApost(index(1)).simData.DataNames,['D_T' TumorArray{k} '_perc']) )];
end

f = figure; hold on; box on;
set(f,'Position', [50 50 800 500]);
% Find a list of input parameter for each PSA case
tempIn = params_in.(param).LHS(index);
paramName = params_in.(param).ScreenName;
paramInModel = sbioselect (model, 'Type', 'parameter', 'Name', param);
if strcmp(paramInModel.ValueUnits,'dimensionless')
    paramUnit = '';
else
    paramUnit = paramInModel.ValueUnits;
end

min_index=NaN(length(index),1);
% Find outputs for PSA cases
j = [find(strcmp(simDataPSApost(index(1)).simData.DataNames,'D_T_all_perc') )];
for i =1:length(index)
    D_T_perc = simDataPSApost(index(i)).simData.Data(:,j);
    % tempOut(i,1) = D_T_perc(end); % end tumor size
    % best overall response assuming tumor measurement every 8 weeks
    [tempOut(i,1),min_index(i,1)] = min(D_T_perc((test_freq+1):test_freq:end));
end

for k=1:ntumors

    for np=1:length(index)

         pos = find(ismember(simDataPSApost(index(np)).simData.DataNames,['D_T' TumorArray{k}]));
         D_T_temp = simDataPSApost(index(np)).simData.Data(:,pos);

         if(D_T_temp(1)>=0.2)
    
            %thein: check timepoint & plot only tumors that exist
            tmp=simDataPSApost(index(np)).simData.Data((test_freq+1):test_freq:end,tumorperc_pos(k,1));
            tumor_response(np,k)=tmp(min_index(np,1));

         end

    end

end

temp = [tempIn, tempOut,tumor_response];
temp = sortrows(temp, 2 ,'descend');
tempMedian = median(  temp(:,1) );

ax = gca;
i = 1:length(index)

ax.ColorOrderIndex = 1; 

%C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]};
C={[0.6350, 0.0780, 0.1840],[0.4660, 0.6740, 0.1880],[0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.3010, 0.7450, 0.9330],[0.4660, 0.6740, 0.1880]};
for k=1:ntumors
    
    hold on;
    h=plot(i,temp(:,2+k),'o','Color',C{k});

end
h2=plot(i,temp(:,2),'-k.');

h3 = plot( [0, n_PSA*1.2], [+20, +20], '--k' );
h4 = plot( [0, n_PSA*1.2], [-30, -30], '--k' );
hx = text(n_PSA*1.10, 60,'PD');
hx = text(n_PSA*1.10,-10,'SD');
hx = text(n_PSA*1.05,-60,'PR/CR');

%legend([h1 h2],{sprintf([paramName,' $<$ %0.2e ',paramUnit,' (median)'],tempMedian),sprintf([paramName,' $>$ %0.2e ',paramUnit,' (median)'],tempMedian)});

% legend([h1 h2],{[param,' < ',num2str(tempMedian,'%6.0f'),' median'],[param,' > ',num2str(tempMedian,'\%6.0f'),' median']});
ylabel('% change in sum of tumor diameters','Fontsize',12);
%ylim([-100 100])
xlim([0  n_PSA*1.2])
set(gca,'Fontsize',14)
set(gca,'XTick',[])
