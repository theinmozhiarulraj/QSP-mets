%% S1 simulations - primary + lung met
%% S2 simulations - only primary
%% S3 simulations - only met

% time points
time=0:2500;

% load the data
s1=load("simDataPSA_s1.mat");
simDataPSA_s1=s1.simDataPSA;

s1=load("params_out_s1.mat");
params_out_s1=s1.params_out;

s2=load("simDataPSA_s2.mat");
simDataPSA_s2=s2.simDataPSA;

s2=load("params_out_s2.mat");
params_out_s2=s2.params_out;

s3=load("simDataPSA_s3.mat");
simDataPSA_s3=s3.simDataPSA;

s3=load("params_out_s3.mat");
params_out_s3=s3.params_out;

% index of simulations that worked in s1, s2 and s3
index_patients=[];
for i=1:length(params_out_s1.patient)

    if(params_out_s1.patient(i)==1 && params_out_s2.patient(i)==1 && params_out_s3.patient(i)==1)
        index_patients=[index_patients,i];
    end

end
disp(index_patients)

npatients = length(index_patients);

% tumor volumes
prim_vol_s1=[];
met_vol_s1=[];

prim_vol_s2=[];
met_vol_s2=[];
 
prim_vol_s3=[];
met_vol_s3=[];


for i = 1:npatients

    index = index_patients(i);
    
    % tumor volumes
    pos = find(ismember(simDataPSA_s1(index).simData.DataNames,'V_T'));
    prim_vol_temp = simDataPSA_s1(index).simData.Data(:,pos);
    prim_vol_s1=[prim_vol_s1,prim_vol_temp];

    prim_vol_temp = simDataPSA_s2(index).simData.Data(:,pos);
    prim_vol_s2=[prim_vol_s2,prim_vol_temp];

    prim_vol_temp = simDataPSA_s3(index).simData.Data(:,pos);
    prim_vol_s3=[prim_vol_s3,prim_vol_temp];

    pos = find(ismember(simDataPSA_s1(index).simData.DataNames,'V_T_Ln1'));
    met_vol_temp = simDataPSA_s1(index).simData.Data(:,pos);
    met_vol_s1=[met_vol_s1,met_vol_temp];

    met_vol_temp = simDataPSA_s2(index).simData.Data(:,pos);
    met_vol_s2=[met_vol_s2,met_vol_temp];

    met_vol_temp = simDataPSA_s3(index).simData.Data(:,pos);
    met_vol_s3=[met_vol_s3,met_vol_temp];


end

%% calculate the difference

difference_prim=prim_vol_s2-prim_vol_s1;
difference_met=met_vol_s3-met_vol_s1;

%% figures
% plot the volume of primary and lung met

% figure 1: primary and met volume in s1
figure();
tiledlayout(10,10);

for j=1:npatients
    nexttile
    plot(time,prim_vol_s1(:,j),'k');
    hold on;
    plot(time,met_vol_s1(:,j),'b');
end

% figure 2: primary and met volume in s2
figure();
tiledlayout(10,10);

for j=1:npatients
    nexttile
    plot(time,prim_vol_s2(:,j),'k');
    hold on;
    plot(time,met_vol_s2(:,j),'b');
end

% figure 3: primary and met volume in s3
figure();
tiledlayout(10,10);

for j=1:npatients
    nexttile
    plot(time,prim_vol_s3(:,j),'k');
    hold on;
    plot(time,met_vol_s3(:,j),'b');

end


% figure 4: tumor volumes in all conditions together

figure();
tiledlayout(10,10);

for j=1:npatients
    nexttile
    plot(time,prim_vol_s1(:,j),'k');
    hold on;
    plot(time,met_vol_s1(:,j),'b');
    hold on;
    plot(time,prim_vol_s2(:,j),':k');
    hold on;
    plot(time,met_vol_s3(:,j),':b');

end

figure();
plot(difference_prim)
hold on;
yline(0,':')
ylim([-10 10])
xlim([0 400])

figure();
plot(difference_met)
hold on;
yline(0,':')
ylim([-10 10])
xlim([0 400])

