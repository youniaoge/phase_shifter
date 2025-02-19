clc;
clear all;

%input
freq=4.4e10;%在哪个频率处进行筛�1ￄ1�7�1�7�1�7�1�7�1�7�1�7
N=6;%%bit敄1�7�1�7�1�7�1�7
j0=65;%参�1ￄ1�7�1�7��1ￄ1�7�1�7�的状�1ￄ1�7�1�7�列号_霄1�7�1�7�1�7�1�7�找Imax+Qmin+的状态点
amp_range_db2=1;%划定的幅度波动的范围(dB)

%f1=figure;
%f2=figure;
f_select_freq=figure;
f_db=figure;
f_phase_absl=figure;
f_phase_rela=figure;
f_rms_phase=figure;
f_rms_db=figure;

%plot const @freq
phase=readtable('PS_30_50_S21_deg_processed.csv'); 
phase_arr=table2array(phase);
k0=find(phase_arr(1:end,1)==freq);
phase_freq_arr=phase_arr(k0,2:end);
amp=readtable('PS_30_50_S21_dB_processed.csv');
amp_arr=table2array(amp);
amp_freq=amp(k0,2:end);
amp_freq_db=table2array(amp_freq);
for n=1:size(amp_freq_db,2)
    amp_freq_mag(1,n)=10^(amp_freq_db(1,n)/20);
end
X_freq=amp_freq_mag.*cos(phase_freq_arr*pi/180);
Y_freq=amp_freq_mag.*sin(phase_freq_arr*pi/180);
%figure(f1);
%scatter(X_freq,Y_freq,8,'blue','filled');

phase_freq_rela=phase_freq_arr-linspace(phase_freq_arr(1,j0),phase_freq_arr(1,j0),size(phase_freq_arr,2));
X_rela=amp_freq_mag.*cos(phase_freq_rela*pi/180);
Y_rela=amp_freq_mag.*sin(phase_freq_rela*pi/180);
%figure(f2);
%scatter(X_rela,Y_rela,8,'blue','filled');


%Select2
mag_max2=amp_freq_mag(1,j0)*(10^(amp_range_db2/20));
mag_min2=amp_freq_mag(1,j0)/(10^(amp_range_db2/20));
X_in2=[];
Y_in2=[];
code_in2=[];
for n=1:size(X_rela,2)
    if mag_min2^2<=(X_rela(1,n)^2+Y_rela(1,n)^2)&&(X_rela(1,n)^2+Y_rela(1,n)^2)<=mag_max2^2
        X_in2=[X_in2,X_rela(1,n)];
        Y_in2=[Y_in2,Y_rela(1,n)];
        code_in2=[code_in2,n];
    end
end

phase_ideal=linspace(-90,360*(2^N-1)/2^N-90,2^N);
X_select2=[];
Y_select2=[];%丄1�7�1�7�1�7�1�7�幅度波动范围内移相精度朄1�7�1�7�1�7�1�7�的
sum_rms_phase=0;
sum_amp=0;
rms_sum_amp=0;
phase_select_tmp = zeros(2^N, 2); % Initialize phase_select array
code_select_tmp=zeros(2^N,1);
for m=1:2^N
    d_phase_min=90;
    for n=1:size(X_in2,2)
        phase_tmp=atan(Y_in2(1,n)/X_in2(1,n))*180/pi;
        if X_in2(1,n)<0
            phase_tmp=phase_tmp+180;
        end
        d_phase_tmp=abs(phase_ideal(1,m)-phase_tmp);
        if d_phase_tmp<d_phase_min
            d_phase_min=d_phase_tmp;
            X_select2(1,m)=X_in2(1,n);
            Y_select2(1,m)=Y_in2(1,n);
            phase_select_tmp(m, :) = [m, phase_tmp]; % Store m and phase_tmp in phase_select
            code_select_tmp(m,1)=code_in2(1,n);
        end
    end
    
    sum_rms_phase=sum_rms_phase+d_phase_min^2;
    sum_amp=sum_amp+sqrt(X_select2(1,m)^2+Y_select2(1,m)^2);
end
rms_phase=sqrt(sum_rms_phase/(2^N-1));
%ave_amp=sum_amp/2^N;
for n=1:size(X_select2,2)
    rms_sum_amp=rms_sum_amp+(sqrt(X_select2(1,n)^2+Y_select2(1,n)^2)-sum_amp/2^N)^2;%这里先拿幅度
end
rms_amp=sqrt(rms_sum_amp/2^N);
figure(f_select_freq)
title('select freq') %在图形的顶端加注文字作为图形各1�7�1�7�1�7�1�7
hold on
viscircles([0,0],mag_max2);
viscircles([0,0],mag_min2);
scatter(X_rela,Y_rela,8,'blue','filled');
scatter(X_in2,Y_in2,12,'green','filled');
scatter(X_select2,Y_select2,16,'red');
hold off

ind_select=[];
for m=1:size(X_select2,2)
    for n=1:size(X_rela,2)
        if X_rela(1,n)==X_select2(1,m)&&Y_rela(1,n)==Y_select2(1,m)
            ind_select=[ind_select,n];
        end
    end
end

%plot phase amp
figure(f_phase_absl)
title('f phase abs') %在图形的顶端加注文字作为图形各1�7�1�7�1�7�1�7
hold on
for m=1:size(ind_select,2)
    plot(phase_arr(1:end,1),phase_arr(1:end,ind_select(1,m)+1));
end
hold off
figure(f_db)
hold on
for m=1:size(ind_select,2)
    plot(phase_arr(1:end,1),amp_arr(1:end,ind_select(1,m)+1));
end
hold off



% figure;
% hold on
% for m=1:size(ind_select,2)
%     plot(phase_arr(1:end,1),phase_arr(1:end,ind_select(1,m)+1)-phase_arr(1:end,j0+1));
% end
% hold off
% [mtmp,jmin]=min(phase_arr(k0,ind_select));
% figure;
% hold on
% for m=1:size(ind_select,2)
%     plot(phase_arr(1:end,1),phase_arr(1:end,ind_select(1,m)+1)-phase_arr(1:end,ind_select(1,jmin)));
% end
% hold off


phase_arr_copy=phase_arr;
for m=1:size(ind_select,2)
    if phase_arr(k0,j0+1)>phase_arr(k0,ind_select(1,m)+1)
        phase_arr_copy(:,ind_select(1,m)+1)=phase_arr_copy(:,ind_select(1,m)+1)+360;
    end
end
figure(f_phase_rela)
hold on
for m=1:size(ind_select,2)
    plot(phase_arr(1:end,1),phase_arr_copy(1:end,ind_select(1,m)+1)-phase_arr(1:end,j0+1));
end
hold off

phase_select=phase_arr_copy(:,ind_select(1,:)+1)-phase_arr(:,j0+1);
[btmp,itmp]=sort(phase_select(k0,:));
phase_select=phase_select(:,itmp);%sort
%plot rms phase error
rms_phase_all=[];
for m=1:size(phase_arr,1)
    rms_phase_tmp=0;
    for n=1:2^N
        rms_phase_tmp=rms_phase_tmp+(phase_select(m,n)-(n-1)*(360/2^N))^2;
    end
    rms_phase_all=[rms_phase_all,sqrt(rms_phase_tmp/(2^N-1))];
end
figure(f_rms_phase)
plot(phase_arr(:,1),rms_phase_all(1,:));
title('rms phase error') %在图形的顶端加注文字作为图形各1�7�1�7�1�7�1�7
rms_db_all=[];
ave_db_all=[];
for m=1:size(amp_arr,1)
    sum_db_tmp=0;
    for n=1:size(ind_select,2)
        sum_db_tmp=sum_db_tmp+amp_arr(m,ind_select(1,n)+1);
    end
    ave_db_all=[ave_db_all,sum_db_tmp/2^N];
end
for m=1:size(amp_arr,1)
    rms_db_tmp=0;
    for n=1:size(ind_select,2)
        rms_db_tmp=rms_db_tmp+(amp_arr(m,ind_select(1,n)+ 1)-ave_db_all(1,m))^2;
    end
    rms_db_all=[rms_db_all,sqrt(rms_db_tmp/2^N)];
end
figure(f_rms_db)
plot(phase_arr(:,1),rms_db_all(1,:));
title('rms gain error') %在图形的顶端加注文字作为图形各1�7�1�7�1�7�1�7


% 将 ind_select 转换为 12 位二进制代码
% Convert to 12-bit binary
binary_codes = dec2bin(code_select_tmp, 12);

% Format binary codes with space between bits 6 and 7
formatted_binary_codes = cell(size(binary_codes, 1), 1);
for i = 1:size(binary_codes, 1)
    formatted_binary_codes{i} = [binary_codes(i, 1:6) '' binary_codes(i, 7:12)];
end

% Write formatted binary codes and phase values to CSV
fileID1 = fopen('combined_array.csv', 'w');
for i = 1:size(formatted_binary_codes, 1)
    fprintf(fileID1, '%s,%f\n', formatted_binary_codes{i}, phase_select_tmp(i, 2));
end
fclose(fileID1);

fileID2 = fopen('output.csv', 'w');
for i = 1:size(binary_codes, 1)
    fprintf(fileID2, '%s\n', binary_codes(i,1:12));
end
fclose(fileID2);



