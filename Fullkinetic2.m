%% =========================================================
%  OPENSIM KINETIC ANALYSIS
%  Inverse Dynamics + Static Optimization
%  Author:ali salehi/M.SC student of sport biomechanics(shiraz university)
% ==========================================================
%% SECTION 1: Impoert data
% =========================================================
%   A:Import Data(Numeric Matrix):
%   ID_data    = double(inverse_dynamics);
%   act_data   = double(LaiUhlrich2022_scaled_StaticOptimization_activation);
%   force_data = double(LaiUhlrich2022_scaled_StaticOptimization_force);
%   B: Import directly from this code
% =========================================================
path_ID    = 'inverse_dynamics.sto';
path_act   = 'LaiUhlrich2022_scaled_StaticOptimization_activation.sto';
path_force = 'LaiUhlrich2022_scaled_StaticOptimization_force.sto';

ID_data    = readmatrix(path_ID,    'FileType','text','NumHeaderLines',7,  'Delimiter','\t');
act_data   = readmatrix(path_act,   'FileType','text','NumHeaderLines',9,  'Delimiter','\t');
force_data = readmatrix(path_force, 'FileType','text','NumHeaderLines',14, 'Delimiter','\t');

% حذف ستون‌های NaN احتمالی
ID_data    = ID_data(:,    all(~isnan(ID_data),    1));
act_data   = act_data(:,   all(~isnan(act_data),   1));
force_data = force_data(:, all(~isnan(force_data), 1));

body_mass = 80;   % ← Athlet's mass
BW        = body_mass * 9.81;

time_id = ID_data(:,1);
fs_id   = round(1/mean(diff(time_id)));

fprintf('✅ Files are loaded :\n');
fprintf('   ID:         %d frames | %d cols\n', size(ID_data,1),    size(ID_data,2));
fprintf('   Activation: %d frames | %d cols\n', size(act_data,1),   size(act_data,2));
fprintf('   Force:      %d frames | %d cols\n', size(force_data,1), size(force_data,2));

%% SECTION 2: Export joint moment from ID
M_hip_r     = ID_data(:,8);
M_hip_l     = ID_data(:,11);
M_hip_add_r = ID_data(:,9);
M_hip_add_l = ID_data(:,12);
M_knee_r    = ID_data(:,17);
M_knee_l    = ID_data(:,19);
M_ankle_r   = ID_data(:,27);
M_ankle_l   = ID_data(:,28);
M_lumbar    = ID_data(:,14);

% Wheight normalisiation
M_hip_r_n   = M_hip_r   / body_mass;
M_hip_l_n   = M_hip_l   / body_mass;
M_knee_r_n  = M_knee_r  / body_mass;
M_knee_l_n  = M_knee_l  / body_mass;
M_ankle_r_n = M_ankle_r / body_mass;
M_ankle_l_n = M_ankle_l / body_mass;
M_lumbar_n  = M_lumbar  / body_mass;

fprintf('✅ Joint Moments استخراج شدند.\n');

%% SECTION 3: Export muscle activation
act_glmax_r   = mean([act_data(:,16), act_data(:,17), act_data(:,18)], 2);
act_glmax_l   = mean([act_data(:,56), act_data(:,57), act_data(:,58)], 2);
act_glmed_r   = mean([act_data(:,19), act_data(:,20), act_data(:,21)], 2);
act_glmed_l   = mean([act_data(:,59), act_data(:,60), act_data(:,61)], 2);
act_vasti_r   = mean([act_data(:,39), act_data(:,40), act_data(:,41)], 2);
act_vasti_l   = mean([act_data(:,79), act_data(:,80), act_data(:,81)], 2);
act_recfem_r  = act_data(:,31);
act_recfem_l  = act_data(:,71);
act_ham_r     = mean([act_data(:,8),  act_data(:,33), act_data(:,34)], 2);
act_ham_l     = mean([act_data(:,48), act_data(:,73), act_data(:,74)], 2);
act_gas_r     = mean([act_data(:,14), act_data(:,15)], 2);
act_gas_l     = mean([act_data(:,54), act_data(:,55)], 2);
act_sol_r     = act_data(:,35);
act_sol_l     = act_data(:,75);
act_tibant_r  = act_data(:,37);
act_tibant_l  = act_data(:,77);
act_hipflex_r = mean([act_data(:,26), act_data(:,30)], 2);
act_hipflex_l = mean([act_data(:,66), act_data(:,70)], 2);

fprintf('✅ Muscle Activations exported.\n');

%% SECTION 4: Export Muscle Forces
force_glmax_r = mean([force_data(:,16), force_data(:,17), force_data(:,18)], 2);
force_glmax_l = mean([force_data(:,56), force_data(:,57), force_data(:,58)], 2);
force_vasti_r = mean([force_data(:,39), force_data(:,40), force_data(:,41)], 2);
force_vasti_l = mean([force_data(:,79), force_data(:,80), force_data(:,81)], 2);
force_gas_r   = mean([force_data(:,14), force_data(:,15)], 2);
force_gas_l   = mean([force_data(:,54), force_data(:,55)], 2);
force_sol_r   = force_data(:,35);
force_sol_l   = force_data(:,75);

fprintf('✅ Muscle Forces exported.\n');

%% SECTION 5: Butterworth filter
[bf,af] = butter(4, 6/(fs_id/2), 'low');

M_hip_r_f     = filtfilt(bf,af,M_hip_r_n);
M_hip_l_f     = filtfilt(bf,af,M_hip_l_n);
M_knee_r_f    = filtfilt(bf,af,M_knee_r_n);
M_knee_l_f    = filtfilt(bf,af,M_knee_l_n);
M_ankle_r_f   = filtfilt(bf,af,M_ankle_r_n);
M_ankle_l_f   = filtfilt(bf,af,M_ankle_l_n);

act_glmax_r_f   = filtfilt(bf,af,act_glmax_r);
act_glmax_l_f   = filtfilt(bf,af,act_glmax_l);
act_glmed_r_f   = filtfilt(bf,af,act_glmed_r);
act_glmed_l_f   = filtfilt(bf,af,act_glmed_l);
act_vasti_r_f   = filtfilt(bf,af,act_vasti_r);
act_vasti_l_f   = filtfilt(bf,af,act_vasti_l);
act_recfem_r_f  = filtfilt(bf,af,act_recfem_r);
act_recfem_l_f  = filtfilt(bf,af,act_recfem_l);
act_ham_r_f     = filtfilt(bf,af,act_ham_r);
act_ham_l_f     = filtfilt(bf,af,act_ham_l);
act_gas_r_f     = filtfilt(bf,af,act_gas_r);
act_gas_l_f     = filtfilt(bf,af,act_gas_l);
act_sol_r_f     = filtfilt(bf,af,act_sol_r);
act_sol_l_f     = filtfilt(bf,af,act_sol_l);
act_tibant_r_f  = filtfilt(bf,af,act_tibant_r);
act_tibant_l_f  = filtfilt(bf,af,act_tibant_l);
act_hipflex_r_f = filtfilt(bf,af,act_hipflex_r);
act_hipflex_l_f = filtfilt(bf,af,act_hipflex_l);

fprintf('✅ .\n');

%% SECTION 6: Distinct Bottom
[~,bottom_id] = max((abs(M_knee_r_f) + abs(M_knee_l_f))/2);

%% SECTION 7: Key vairiables
fprintf('\n=== Joint Moments (N.m/kg) ===\n');
fprintf('%-35s %-10s %-10s\n','Variable','Right','Left');
fprintf('%s\n',repmat('-',1,57));
fprintf('%-35s %-10.2f %-10.2f\n','Peak Hip Moment',   max(abs(M_hip_r_f)),   max(abs(M_hip_l_f)));
fprintf('%-35s %-10.2f %-10.2f\n','Peak Knee Moment',  max(abs(M_knee_r_f)),  max(abs(M_knee_l_f)));
fprintf('%-35s %-10.2f %-10.2f\n','Peak Ankle Moment', max(abs(M_ankle_r_f)), max(abs(M_ankle_l_f)));
fprintf('%-35s %-10.2f\n',        'Peak Lumbar Moment',max(abs(M_lumbar_n)));

fprintf('\n=== Peak Muscle Activations ===\n');
fprintf('%-25s %-10s %-10s\n','Muscle','Right','Left');
fprintf('%s\n',repmat('-',1,47));
fprintf('%-25s %-10.3f %-10.3f\n','Gluteus Maximus',   max(act_glmax_r_f),   max(act_glmax_l_f));
fprintf('%-25s %-10.3f %-10.3f\n','Gluteus Medius',    max(act_glmed_r_f),   max(act_glmed_l_f));
fprintf('%-25s %-10.3f %-10.3f\n','Vasti',             max(act_vasti_r_f),   max(act_vasti_l_f));
fprintf('%-25s %-10.3f %-10.3f\n','Rectus Femoris',    max(act_recfem_r_f),  max(act_recfem_l_f));
fprintf('%-25s %-10.3f %-10.3f\n','Hamstrings',        max(act_ham_r_f),     max(act_ham_l_f));
fprintf('%-25s %-10.3f %-10.3f\n','Gastrocnemius',     max(act_gas_r_f),     max(act_gas_l_f));
fprintf('%-25s %-10.3f %-10.3f\n','Soleus',            max(act_sol_r_f),     max(act_sol_l_f));
fprintf('%-25s %-10.3f %-10.3f\n','Tibialis Anterior', max(act_tibant_r_f),  max(act_tibant_l_f));
fprintf('%-25s %-10.3f %-10.3f\n','Hip Flexors',       max(act_hipflex_r_f), max(act_hipflex_l_f));

SI_hip_m   = abs(max(abs(M_hip_r_f))  -max(abs(M_hip_l_f)))  /((max(abs(M_hip_r_f))  +max(abs(M_hip_l_f)))/2)  *100;
SI_knee_m  = abs(max(abs(M_knee_r_f)) -max(abs(M_knee_l_f))) /((max(abs(M_knee_r_f)) +max(abs(M_knee_l_f)))/2) *100;
SI_ankle_m = abs(max(abs(M_ankle_r_f))-max(abs(M_ankle_l_f)))/((max(abs(M_ankle_r_f))+max(abs(M_ankle_l_f)))/2)*100;

fprintf('\n=== Moment Symmetry Index ===\n');
fprintf('Hip:   %.1f%%\n',SI_hip_m);
fprintf('Knee:  %.1f%%\n',SI_knee_m);
fprintf('Ankle: %.1f%%\n',SI_ankle_m);

%% SECTION 8: Plot 1 — Joint Moments
figure('Name','Joint Moments','Position',[50 50 1400 600]);

subplot(1,3,1);
plot(time_id,M_hip_r_f,'b-',time_id,M_hip_l_f,'r--','LineWidth',1.5);
yline(0,'k--'); xline(time_id(bottom_id),'k:');
xlabel('Time (s)'); ylabel('Moment (N·m/kg)');
title('Hip Extension Moment'); legend('Right','Left'); grid on;
subplot(1,3,2);
plot(time_id,M_knee_r_f,'b-',time_id,M_knee_l_f,'r--','LineWidth',1.5);
yline(0,'k--'); xline(time_id(bottom_id),'k:');
xlabel('Time (s)'); ylabel('Moment (N·m/kg)');
title('Knee Extension Moment'); legend('Right','Left'); grid on;

subplot(1,3,3);
plot(time_id,M_ankle_r_f,'b-',time_id,M_ankle_l_f,'r--','LineWidth',1.5);
yline(0,'k--'); xline(time_id(bottom_id),'k:');
xlabel('Time (s)'); ylabel('Moment (N·m/kg)');
title('Ankle Moment'); legend('Right','Left'); grid on;

sgtitle('Net Joint Moments — OpenSim Inverse Dynamics','FontSize',13,'FontWeight','bold');

%% SECTION 9: Plot 2 — Muscle Activations
figure('Name','Muscle Activations','Position',[50 50 1400 900]);

muscles_R  = {act_glmax_r_f,act_glmed_r_f,act_vasti_r_f,act_recfem_r_f, ...
              act_ham_r_f,act_gas_r_f,act_sol_r_f,act_tibant_r_f,act_hipflex_r_f};
muscles_L  = {act_glmax_l_f,act_glmed_l_f,act_vasti_l_f,act_recfem_l_f, ...
              act_ham_l_f,act_gas_l_f,act_sol_l_f,act_tibant_l_f,act_hipflex_l_f};
mus_names  = {'Glut Max','Glut Med','Vasti','Rect Fem', ...
              'Hamstrings','Gastrocnemius','Soleus','Tib Ant','Hip Flexors'};

for i = 1:9
    subplot(3,3,i);
    plot(time_id,muscles_R{i},'b-','LineWidth',1.5); hold on;
    plot(time_id,muscles_L{i},'r--','LineWidth',1.5);
    xline(time_id(bottom_id),'k:');
    ylim([0 1.05]); xlabel('Time (s)'); ylabel('Activation');
    title(mus_names{i}); legend('Right','Left','FontSize',7); grid on;
end
sgtitle('Muscle Activations — Static Optimization','FontSize',13,'FontWeight','bold');

%% SECTION 10: Plot 3 — Activation Heatmap
figure('Name','Heatmap','Position',[50 50 1200 600]);

act_mat_R = cell2mat(cellfun(@(x) x', muscles_R, 'UniformOutput', false));
act_mat_L = cell2mat(cellfun(@(x) x', muscles_L, 'UniformOutput', false));

subplot(1,2,1);
imagesc(time_id, 1:9, act_mat_R); colormap(hot); colorbar;
set(gca,'YTickLabel',mus_names,'YTick',1:9);
xlabel('Time (s)'); title('Activation Heatmap — Right'); clim([0 1]);

subplot(1,2,2);
imagesc(time_id, 1:9, act_mat_L); colormap(hot); colorbar;
set(gca,'YTickLabel',mus_names,'YTick',1:9);
xlabel('Time (s)'); title('Activation Heatmap — Left'); clim([0 1]);

sgtitle('Muscle Activation Heatmap','FontSize',13,'FontWeight','bold');

%% SECTION 11: Plot 4 — Muscle Forces
figure('Name','Muscle Forces','Position',[50 50 1200 600]);

subplot(2,2,1);
plot(time_id,force_glmax_r,'b-',time_id,force_glmax_l,'r--','LineWidth',1.5);
xline(time_id(bottom_id),'k:');
xlabel('Time (s)'); ylabel('Force (N)');
title('Gluteus Maximus'); legend('Right','Left'); grid on;

subplot(2,2,2);
plot(time_id,force_vasti_r,'b-',time_id,force_vasti_l,'r--','LineWidth',1.5);
xline(time_id(bottom_id),'k:');
xlabel('Time (s)'); ylabel('Force (N)');
title('Vasti'); legend('Right','Left'); grid on;

subplot(2,2,3);
plot(time_id,force_gas_r,'b-',time_id,force_gas_l,'r--','LineWidth',1.5);
xline(time_id(bottom_id),'k:');
xlabel('Time (s)'); ylabel('Force (N)');
title('Gastrocnemius'); legend('Right','Left'); grid on;

subplot(2,2,4);
plot(time_id,force_sol_r,'b-',time_id,force_sol_l,'r--','LineWidth',1.5);
xline(time_id(bottom_id),'k:');
xlabel('Time (s)'); ylabel('Force (N)');
title('Soleus'); legend('Right','Left'); grid on;

sgtitle('Muscle Forces — Static Optimization','FontSize',13,'FontWeight','bold');

%% SECTION 12: Plot 5 — Symmetry
figure('Name','Kinetic Symmetry','Position',[50 50 900 500]);

subplot(1,2,1);
SI_m = [SI_hip_m, SI_knee_m, SI_ankle_m];
b_si = bar(SI_m,'FaceColor','flat');
for i=1:3
    if SI_m(i)<=10;     b_si.CData(i,:)=[0.2 0.7 0.3];
    elseif SI_m(i)<=20; b_si.CData(i,:)=[1.0 0.6 0.0];
    else;               b_si.CData(i,:)=[0.9 0.2 0.2];
    end
end
yline(10,'k--','10%');
set(gca,'XTickLabel',{'Hip','Knee','Ankle'});
ylabel('SI%'); title('Moment Symmetry Index'); grid on;
subplot(1,2,2);
SI_act = zeros(1,9);
for i=1:9
    R=max(muscles_R{i}); L=max(muscles_L{i});
    if (R+L)>1e-6; SI_act(i)=abs(R-L)/((R+L)/2)*100; end
end
b_act = bar(SI_act,'FaceColor','flat');
for i=1:9
    if SI_act(i)<=10;     b_act.CData(i,:)=[0.2 0.7 0.3];
    elseif SI_act(i)<=20; b_act.CData(i,:)=[1.0 0.6 0.0];
    else;                 b_act.CData(i,:)=[0.9 0.2 0.2];
    end
end
yline(10,'k--','10%');
set(gca,'XTickLabel',mus_names,'XTickLabelRotation',30,'FontSize',7);
ylabel('SI%'); title('Activation Symmetry Index'); grid on;

sgtitle('Kinetic Symmetry Analysis','FontSize',13,'FontWeight','bold');

%% SECTION 13: Save
T = table({'Hip';'Knee';'Ankle'}, ...
    [max(abs(M_hip_r_f)); max(abs(M_knee_r_f)); max(abs(M_ankle_r_f))], ...
    [max(abs(M_hip_l_f)); max(abs(M_knee_l_f)); max(abs(M_ankle_l_f))], ...
    [SI_hip_m; SI_knee_m; SI_ankle_m], ...
    'VariableNames',{'Joint','Right_Nm_kg','Left_Nm_kg','SI_percent'});
writetable(T,'opensim_kinetics_results.csv');
fprintf('\n✅ نتایج ذخیره شد.\n');
fprintf('⚠️  Activations خیلی بالا — با احتیاط interpret کن.\n');