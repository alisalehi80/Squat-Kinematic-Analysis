%% =========================================================
%  This pipeline imports,filters,and processes full lower-limb kinematic
%  data
%  Author: ali_salehi/M.SC student of sport Biomechanics(shirazuniversity)
%  Date:June 2026
%  This code need Opencap kinematic data(mot) as numeric matrixe
%  SQUAT KINEMATIC ANALYSIS
%  Data source: Squat.mot | Header Lines: 11 | Numeric Matrix
% ==========================================================

clc; close all;

%% SECTION 1: import data
raw  = double(Squat);
data = raw;

for col = 1:size(data,2)
    col_data = data(:,col);
    nan_idx  = isnan(col_data) | ~isfinite(col_data);
    if any(nan_idx)
        t_all  = (1:length(col_data))';
        t_good = t_all(~nan_idx);
        v_good = col_data(~nan_idx);
        if length(t_good) > 1
            data(:,col) = interp1(t_good,v_good,t_all,'linear','extrap');
        end
    end
end

if any(~isfinite(data(:)))
    warning('we have Nan still');
else
    fprintf('✅ داده‌ها finite هستند.\n');
end

time = data(:,1);
fs   = 60;
dt   = 1/fs;
n    = size(data,1);
fprintf(' داده بارگذاری شد | فریم: %d | زمان: %.2f s\n', n, time(end));

%% SECTION 2: Extracting angles
hip_flex_R  = data(:,8);
hip_flex_L  = data(:,15);
hip_add_R   = data(:,9);
hip_add_L   = data(:,16);
knee_flex_R = data(:,11);
knee_flex_L = data(:,18);
ankle_df_R  = data(:,12);
ankle_df_L  = data(:,19);
lumbar_ext  = data(:,22);
pelvis_tilt = data(:,2);
pelvis_list = data(:,3);
fprintf('✅ Angles extracted.\n');

%% SECTION 3: Filtering
fc = 6;
[b,a] = butter(4, fc/(fs/2), 'low');

hip_flex_R_f  = filtfilt(b,a,hip_flex_R);
hip_flex_L_f  = filtfilt(b,a,hip_flex_L);
hip_add_R_f   = filtfilt(b,a,hip_add_R);
hip_add_L_f   = filtfilt(b,a,hip_add_L);
knee_flex_R_f = filtfilt(b,a,knee_flex_R);
knee_flex_L_f = filtfilt(b,a,knee_flex_L);
ankle_df_R_f  = filtfilt(b,a,ankle_df_R);
ankle_df_L_f  = filtfilt(b,a,ankle_df_L);
lumbar_ext_f  = filtfilt(b,a,lumbar_ext);
pelvis_tilt_f = filtfilt(b,a,pelvis_tilt);
pelvis_list_f = filtfilt(b,a,pelvis_list);
fprintf('✅ فیلتر اعمال شد.\n');

%% SECTION 4: Squat phases
knee_avg       = (knee_flex_R_f + knee_flex_L_f)/2;
[~,bottom_idx] = max(knee_avg);
fprintf('\n=== Squatphases ===\n');
fprintf('Bottom: Frame %d | زمان %.3f s\n', bottom_idx, time(bottom_idx));
fprintf('پایین رفتن: %.3f s | بالا اومدن: %.3f s\n', ...
    time(bottom_idx)-time(1), time(end)-time(bottom_idx));

%% SECTION 5: Velocity و Acceleration
hip_vel_R   = gradient(hip_flex_R_f,  dt);
hip_vel_L   = gradient(hip_flex_L_f,  dt);
knee_vel_R  = gradient(knee_flex_R_f, dt);
knee_vel_L  = gradient(knee_flex_L_f, dt);
ankle_vel_R = gradient(ankle_df_R_f,  dt);
ankle_vel_L = gradient(ankle_df_L_f,  dt);
hip_acc_R   = gradient(hip_vel_R,  dt);
hip_acc_L   = gradient(hip_vel_L,  dt);
knee_acc_R  = gradient(knee_vel_R, dt);
knee_acc_L  = gradient(knee_vel_L, dt);
ankle_acc_R = gradient(ankle_vel_R,dt);
ankle_acc_L = gradient(ankle_vel_L,dt);

%% SECTION 6: Variables
hip_ROM_R   = range(hip_flex_R_f);
hip_ROM_L   = range(hip_flex_L_f);
knee_ROM_R  = range(knee_flex_R_f);
knee_ROM_L  = range(knee_flex_L_f);
ankle_ROM_R = range(ankle_df_R_f);
ankle_ROM_L = range(ankle_df_L_f);

fprintf('\n=== Kinematic variables ===\n');
fprintf('%-40s %-10s %-10s\n','Variable','Right','Left');
fprintf('%s\n',repmat('-',1,62));
fprintf('%-40s %-10.1f %-10.1f\n','Hip Flexion ROM (deg)',        hip_ROM_R,                 hip_ROM_L);
fprintf('%-40s %-10.1f %-10.1f\n','Knee Flexion ROM (deg)',       knee_ROM_R,                knee_ROM_L);
fprintf('%-40s %-10.1f %-10.1f\n','Ankle Dorsiflexion ROM (deg)', ankle_ROM_R,               ankle_ROM_L);
fprintf('%-40s %-10.1f %-10.1f\n','Peak Hip Flex at Bottom',      hip_flex_R_f(bottom_idx),  hip_flex_L_f(bottom_idx));
fprintf('%-40s %-10.1f %-10.1f\n','Peak Knee Flex at Bottom',     knee_flex_R_f(bottom_idx), knee_flex_L_f(bottom_idx));
fprintf('%-40s %-10.1f %-10.1f\n','Peak Ankle DorsiFlex',         ankle_df_R_f(bottom_idx),  ankle_df_L_f(bottom_idx));
fprintf('%-40s %-10.1f %-10.1f\n','Peak Hip Adduction (deg)',     max(abs(hip_add_R_f)),     max(abs(hip_add_L_f)));fprintf('%-40s %-10.1f\n',        'Peak Lumbar Extension (deg)',  max(abs(lumbar_ext_f)));
fprintf('%-40s %-10.1f\n',        'Peak Pelvis Tilt (deg)',       max(abs(pelvis_tilt_f)));
fprintf('%-40s %-10.1f\n',        'Peak Pelvis List (deg)',       max(abs(pelvis_list_f)));
fprintf('%-40s %-10.1f %-10.1f\n','Peak Knee Ang Vel (deg/s)',    max(abs(knee_vel_R)),       max(abs(knee_vel_L)));
fprintf('%-40s %-10.1f %-10.1f\n','Peak Hip Ang Vel (deg/s)',     max(abs(hip_vel_R)),        max(abs(hip_vel_L)));
fprintf('%-40s %-10.1f %-10.1f\n','Peak Ankle Ang Vel (deg/s)',   max(abs(ankle_vel_R)),      max(abs(ankle_vel_L)));
fprintf('%-40s %-10.1f %-10.1f\n','Peak Knee Ang Acc (deg/s²)',   max(abs(knee_acc_R)),       max(abs(knee_acc_L)));
fprintf('%-40s %-10.1f %-10.1f\n','Peak Hip Ang Acc (deg/s²)',    max(abs(hip_acc_R)),        max(abs(hip_acc_L)));
fprintf('%-40s %-10.1f %-10.1f\n','Peak Ankle Ang Acc (deg/s²)',  max(abs(ankle_acc_R)),      max(abs(ankle_acc_L)));

%% SECTION 7: Symmetry Index
var_names = {'Hip ROM','Knee ROM','Ankle ROM','Hip Peak','Knee Peak','Ankle Peak'};
vars_R    = [hip_ROM_R, knee_ROM_R, ankle_ROM_R, ...
             hip_flex_R_f(bottom_idx), knee_flex_R_f(bottom_idx), ankle_df_R_f(bottom_idx)];
vars_L    = [hip_ROM_L, knee_ROM_L, ankle_ROM_L, ...
             hip_flex_L_f(bottom_idx), knee_flex_L_f(bottom_idx), ankle_df_L_f(bottom_idx)];
SI_values = zeros(1,6);
fprintf('\n=== Symmetry Index ===\n');
fprintf('%-30s %-10s %s\n','Variable','SI%','Interpretation');
fprintf('%s\n',repmat('-',1,55));
for i = 1:6
    denom = (abs(vars_R(i))+abs(vars_L(i)))/2;
    if denom < 1e-6; SI_values(i) = 0;
    else; SI_values(i) = abs(vars_R(i)-vars_L(i))/denom*100;
    end
    if SI_values(i)<=10;      s='✅ Symmetric';
    elseif SI_values(i)<=20;  s='⚠️  Mild Asymmetry';
    else;                     s='❌ Significant Asymmetry';
    end
    fprintf('%-30s %-10.1f %s\n',var_names{i},SI_values(i),s);
end

%% SECTION 8: Plots
figure('Name','Squat Kinematic Analysis','Position',[50 50 1400 900]);
subplot(3,3,1);
plot(time,hip_flex_R_f,'b-',time,hip_flex_L_f,'r--','LineWidth',1.5);
xline(time(bottom_idx),'k--'); xlabel('Time (s)'); ylabel('deg');
title('Hip Flexion'); legend('Right','Left'); grid on;
subplot(3,3,2);
plot(time,knee_flex_R_f,'b-',time,knee_flex_L_f,'r--','LineWidth',1.5);
xline(time(bottom_idx),'k--'); xlabel('Time (s)'); ylabel('deg');
title('Knee Flexion'); legend('Right','Left'); grid on;
subplot(3,3,3);
plot(time,ankle_df_R_f,'b-',time,ankle_df_L_f,'r--','LineWidth',1.5);
xline(time(bottom_idx),'k--'); xlabel('Time (s)'); ylabel('deg');
title('Ankle Dorsiflexion'); legend('Right','Left'); grid on;
subplot(3,3,4);
plot(time,hip_vel_R,'b-',time,hip_vel_L,'r--','LineWidth',1.5);
xlabel('Time (s)'); ylabel('deg/s'); title('Hip Angular Velocity'); grid on;
subplot(3,3,5);
plot(time,knee_vel_R,'b-',time,knee_vel_L,'r--','LineWidth',1.5);
xlabel('Time (s)'); ylabel('deg/s'); title('Knee Angular Velocity'); grid on;
subplot(3,3,6);
plot(time,ankle_vel_R,'b-',time,ankle_vel_L,'r--','LineWidth',1.5);
xlabel('Time (s)'); ylabel('deg/s'); title('Ankle Angular Velocity'); grid on;
subplot(3,3,7);
plot(time,lumbar_ext_f,'k-','LineWidth',1.5); hold on;
plot(time,pelvis_tilt_f,'m--','LineWidth',1.5);
plot(time,pelvis_list_f,'g:','LineWidth',1.5);
xline(time(bottom_idx),'k--');
xlabel('Time (s)'); ylabel('deg');
title('Trunk & Pelvis');
legend('Lumbar Ext','Pelvis Tilt','Pelvis List'); grid on;
subplot(3,3,8);
plot(time,hip_add_R_f,'b-',time,hip_add_L_f,'r--','LineWidth',1.5);
xline(time(bottom_idx),'k--');
xlabel('Time (s)'); ylabel('deg');
title('Hip Adduction (Valgus)'); legend('Right','Left'); grid on;
subplot(3,3,9);
b_h = bar(SI_values,'FaceColor','flat');
for i=1:6
    if SI_values(i)<=10;     b_h.CData(i,:)=[0.2 0.7 0.3];elseif SI_values(i)<=20; b_h.CData(i,:)=[1.0 0.6 0.0];
    else;                    b_h.CData(i,:)=[0.9 0.2 0.2];
    end
end
yline(10,'k--');
set(gca,'XTickLabel',var_names,'XTickLabelRotation',30,'FontSize',8);
ylabel('SI%'); title('Symmetry Index'); grid on;
sgtitle('Squat Kinematic Analysis — OpenSim IK','FontSize',13,'FontWeight','bold');

writetable(table(var_names',vars_R',vars_L',SI_values', ...
    'VariableNames',{'Variable','Right_deg','Left_deg','SI_percent'}), ...
    'squat_kinematics_results.csv');