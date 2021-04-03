clear all
% Load file am thanh voi bien do ghi vao mang a
% Tan so lay mau Fs
% [x,Fs] = audioread('lab_male.wav');
% [x,Fs] = audioread('lab_female.wav');
% [x,Fs] = audioread('studio_male.wav');
 [x,Fs] = audioread('studio_female.wav');

%===========================CHUONG TRINH CHINH=============================
% Do dai moi khung tin hieu
fr_time = 0.02;
% So luong mau tren moi khung tin hieu
fr_len = fr_time*Fs;
% Tong so luong khung hinh
fr_num = floor(length(x)/fr_len);


%--------------------------Thuat toan Tu tuong quan------------------------

% Chuyen gia tri x ve dang ma tran fr_num hang, fr_len cot 
data_fr = convertx(fr_num,fr_len,x);
% Ham tinh gia tri cua so hamming
w = hamming_window(fr_len);
% Ham tim F0
[F0_TTQ, xcorr_fr, des_maxpeak_fr, delay_fr]  = Find_F0(fr_num,fr_len,x,Fs,w);
% Loc trung vi voi N = 7;
 F0_TTQ = medfilt1(F0_TTQ,7);
% Ham tinh tan so co ban trung binh cua file
F0_ave_TTQ = funct_F0_ave(F0_TTQ);

figure(1);
subplot(211);
plot(x);
title('Tin hieu dau vao','LineWidth',1);
xlabel('Time (s)');
ylabel('Amplitude ');

subplot(212);
t = 1:fr_len:fr_len*fr_num;
stem(t,F0_TTQ,'.', 'r','linestyle','none');
title('Duong bieu dien tan so co ban F0, Ham tu tuong quan');
xlabel('Delay');
ylabel('F0(Hz) ');

figure(2);
des = 50;
subplot(211);
plot(data_fr(des,:));
title('Do thi khung tin hieu do dai 20ms','LineWidth',1);
xlabel('Time(s)');
ylabel('Amplitude');

if( F0_TTQ(des) == 0) 
    title1 = 'Do thi cua khung tin hieu khong tuan hoan';
else 
    title1 =['Do thi ham tu tuong quan cua tin hieu tuan hoan co F0 = ', num2str(F0_TTQ(des))];
end
subplot(212);
plot(xcorr_fr(des,:));
yline(0.3*max(xcorr_fr(des,:)),"Color","red");
title(title1,'LineWidth',1);
hold on;
if(F0_TTQ(des)~= 0)
    plot(delay_fr(des), des_maxpeak_fr(des),"v");
end 
xlabel('Delay');
ylabel('Amplitude');
hold off;


%--------------------------------------------------------------------------

%----------------------------Thuat toan vi sai bien do---------------------

% Sao chep data x sang data u de dung thuat toan amdf
u = repelem(x,1);
% Xu ly tin hieu u de thuc hien thuat toan amdf
u =sigprocessing(x);
% Chuyen gia tri x ve dang ma tran fr_num hang, fr_len cot
data_fr = convertx(fr_num,fr_len,u);
% Ham tim F0
[F0_AMDF, amdf_fr] = Find_F0_amdf(fr_num,fr_len,u,Fs);
% Loc trung vi voi N = 7;
F0_AMDF = medfilt1(F0_AMDF,7);
% Ham tinh gia tri trung binh tan so co ban F0 
F0_ave_AMDF = funct_F0_ave(F0_AMDF);

figure(3);
subplot(211);
plot(x);
subplot(212);
t = 1:fr_len:fr_len*fr_num;
% F0_AMDF = medfilt1(F0_AMDF,5);
stem(t,F0_AMDF,'.', 'r','linestyle','none');
title('Duong bieu dien tan so co ban F0, Ham AMDF');
xlabel('Delay');
ylabel('F0(Hz) ');

figure(4);
des = 120;
subplot(211);
plot(data_fr(des,:));
title('Do thi khung tin hieu do dai 20ms','LineWidth',1);
xlabel('Time(s)');
ylabel('Amplitude');

subplot(212);
plot(amdf_fr(des,:));
title('Do thi khung tin hieu tuan hoan, ham AMDF','LineWidth',1);
hold on;
xlabel('Delay');
ylabel('Amplitude');
hold off;

figure (5);
des1 = 65;
subplot(211);
plot(data_fr(des1,:));
title('Do thi khung tin hieu do dai 20ms','LineWidth',1);
xlabel('Time(s)');
ylabel('Amplitude');

subplot(212);
plot(amdf_fr(des1,:));
title('Do thi khung tin hieu khong tuan hoan, ham AMDF','LineWidth',1);
hold on;
xlabel('Delay');
ylabel('Amplitude');
hold off;
%--------------------------------------------------------------------------

%----------------------------Thuat toan bien doi FFT-----------------------
figure(6);
% Xuat tin hieu ra o mien thoi gian 
t=1/Fs:1/Fs:(length(x)/Fs);     % Truc thoi gian
subplot(211), plot(t,x,'LineWidth',1);
title('Tin hieu ban dau');
xlabel('Thoi gian (s)');
ylabel('Bien do');

fr_time = 0.02; %fr_time(s): Do dai moi khung hinh

fr_len = fr_time*Fs; %fr_len: So mau tren moi khung hinh

N = length(x); %N: tong so mau cua tin hieu

% Tinh so luong khung cua tin hieu (30ms moi khung)
fr_num = floor(N/fr_len);

% Nhan tin hieu voi cua so Hamming
x = x.*hamming(N);

% Tao vecto khung chua tin hieu co do dai 20ms
khung = zeros(1,fr_len);

N_FFT = 16192;% So luong mau tan so(luy thua cua 2)

% Roi rac hoa truc tan so
k=1:N_FFT;
w=k*Fs/N_FFT;

% Tim khung co nang luong lon nhat, thuc hien tinh toan tren khung do
% PowMax Nang luong cuc dai, Start vi tri bat dau cua khung co PowMax
PowMax = 0; Start = 1;
for i = 1:fr_num
    
    for j =(i-1)*fr_len +1 : (fr_len*i)
        khung(j) = x(j);        
    end  
    
    Pow = 0;    % Bien tam thoi, luu gia tri nang luong cua khung dang duoc xet den  
    for j = (i-1)*fr_len +1 : (fr_len*i)
        Pow = Pow + khung(j).*2;
    end    
    
    % Tim khung chua nang luong lon nhat cung diem bat dau cua no
    if Pow > PowMax
        PowMax = Pow;
        Start = i;
    end
end

% Signal la tin hieu co nang luong cuc dai, chinh la tieng noi
% Su dung Signal de xac dinh tan so co ban F0
Signal = x((Start-1)*fr_len +1 : (fr_len*Start));

% Bien doi Fast Fourier Transform tren Signal, voi 
FFT = abs(fft(Signal, N_FFT));


% Mang magnitude chua gia tri bien do cua Signal da duoc bien doi FFT
% Mang co do dai tu 1-350, moi chi so tuong ung voi gia tri cua tan so
magnitude = (FFT(1:350));

% Gia tri ban dau cua tan so
F0 = 75; 

% Gia tri bien do t?i tan so bat dau
max = magnitude(110);

% Tai diem co gia tri bien do cuc dai,
% Xem do chinh la tan so co ban F0 can tim
for i = 110:length(magnitude)
    if max < magnitude(i)
         max = magnitude(i);
         F0 = i;
    end
end

% Ve FFT ra o mien tan so
subplot(212), plot(w(1:N_FFT/2),FFT(1:N_FFT/2));
xlabel('Tan so (Hz)');
ylabel('Bien do');
hold on;
% Ve ra gioi han gia tri tan so tu 75Hz den 350Hz
xline((75), '--r',1), xline((350), '--r',1);
% Danh dau dinh pho ma ta xac nhan do la tan so F0 
hold on;
plot(w(F0), FFT(F0), 'x');

% Gia tri cua thuc cua tan so
F0 = F0*Fs/N_FFT;
title(['Khung tin hieu co nang luong cuc dai, Tan so F0 = ',num2str(F0),' (Hz)']);
% In ket qua F0
fprintf('Tan so co ban F0: %0.2f Hz\n', F0);
%--------------------------------------------------------------------------

%==========================DINH NGHIA HAM==================================

% Thuat toan dua du lieu x ve dang ma tran
% Dau vao: gia tri cua tin hieu x
% Dau ra: Ma tran chua gia tri tren tung khung hinh
function [data_fr] = convertx(fr_num,fr_len,x)
    for i = 0:fr_num-2
        temp = x(i*fr_len+1:(i+1)*fr_len);
        data_fr(i+1,:) = temp;
    end 
 end 

% Ham tim tan so co ban cua tung khung va diem co bien do cuc dai
% fr_num: tong so khung cua tin hieu
% fr_len: tong so tin hieu tren moi khung
% x: tin hieu dau vao sau khi loc trung vi
% Fs: tan so lay mau cua tin hieu
function [F0, xcorr_fr, des_maxpeak_fr, delay_fr] = Find_F0(fr_num,fr_len,x,Fs,w)
    xcorr_fr = zeros(fr_num,fr_len);
    for k = 1 : fr_num
        fr_index = (k-1)*fr_len; %chi so mang cua tin hieu x tai khung thu t
        r = zeros(1,fr_len);
    for delay = 1 : fr_len 
        for j = 1 : fr_len - delay
            r(delay) = r(delay) + x(fr_index + j) * x(fr_index + j + delay - 1)* w(j)^2;
            xcorr_fr(k,delay) = xcorr_fr(k,delay)+ x(fr_index + j) * x(fr_index + j + delay -1 )* w(j)^2;
        end                                     % Cong thuc ham tu tuong quan
    end
         if (r(1) < 0.03)                       % Loai bo nhung khung co gia tri dau tien be hon nguong 0.03
            F = 0;
         else
            F = 0;
            max_peak = 0;
            delay = 0;
            for i = 2 : length(r) - 1
                if (r(i) > r(i -1) && r(i) > r(i+1) && r(i) >= max_peak)   % Ham tim diem cuc dai
                    max_peak = r(i);                % Lay diem cuc dai 
                    delay = i;
                end
            end
            if ( max_peak >= xcorr_fr(1,1)*0.3) % Lay nguong 30% diem cuc dai
                F = Fs/delay;                   % Tinh F0
                if (F <= 80 || F >= 400)        % So sanh dieu kien thoa man F0
                    F = 0;
                end
            end
            delay_fr(k) = delay;                % Luu vi tri delay
            des_maxpeak_fr(k) = max_peak;       % Luu vi tri cac dinh delay
        end
        F0(k) = F;
    end
end
% Gia tri cua so hamming
% Dau vao: do dai khung cua so (length_w)
%-Dau ra: gia tri ham cua so hamming
 function [w]= hamming_window(fr_len)
     w=[];
     for i=1:fr_len
         w(i) = 0.54 - 0.46*cos(2*pi*(i-1)/(fr_len-1));    % cong thuc hamming
     end
 end
    

% Ham tim tan so co ban cua tung khung va diem co bien do cuc dai
% fr_num: tong so khung cua tin hieu
% fr_len: tong so tin hieu tren moi khung
% x: tin hieu dau vao sau khi loc trung vi
% Fs: tan so lay mau cua tin hieu
function [F0, amdf_fr] = Find_F0_amdf(fr_num,fr_len,x,Fs)
    amdf_fr = zeros(fr_num,fr_len);
    for k = 1 : fr_num
        fr_index = (k-1)*fr_len; %chi so mang cua tin hieu x tai khung thu k
        for delay = 1 : fr_len 
            for j = 1 : fr_len - delay
                amdf_fr(k,delay) = amdf_fr(k,delay)+ abs(x(fr_index + j) - x(fr_index + j + delay -1 ));
            end
        end  
       [pks, y] = findpeaks(-amdf_fr(k,40:190));  % Tim cac dinh sau khi tinh amdf tren 1 khung     
        max_peak = max(pks);                       % Tim dinh lon nhat 
        F = 0;
        delay = 0;
        for i = 1: length(pks)
            if(max_peak == pks(i))           % So sanh dinh lon nhat voi mang dinh de lay vi tri tre
                delay = y(i);                 % vi tri tre
                F = Fs/delay;                  % Tim ta so co ban
            end
        end 
        if ( F <= 80 || F >= 400)               % kiem tra tan so co ban thuoc khoang 80-400Hz
            F = 0;                              % Gan vao mang gia tri F0
        end 
        F0(k) = F;         
    end 
end

% Ham xu ly tin hieu data de dung cho thuat toan amdf
% Dau vao: Mang data x
% Dau ra: Mang data u sau khi xu ly
function [u] =sigprocessing(x)
for i =1:length(x)
    u(i) = abs(x(i));       % Lay tri gia tri cua tin hieu tai mau thu i
    if(u(i) < 0.004)        % Neu mau thu i < 0.004 thi cho gia tri mau = 0
        u(i) = 0;
    end 
end 
end 
function [amdf_fr] = AMDF(fr_num,fr_len,x)
    amdf_fr = zeros(fr_num,fr_len);
    for k = 1 : fr_num
        fr_index = (k-1)*fr_len;
        for delay = 1 : fr_len 
            for j = 1 : fr_len - delay
                amdf_fr(k,delay) = amdf_fr(k,delay)+ abs(x(fr_index + j) - x(fr_index + j + delay -1 ));
            end
        end
    end 
end 

% Ham tim gia tri trung binh cua F0
% Dau vao: Mang gia tri cua F0
% Dau ra: Gia tri trung binh F0 cua file tin hieu
function [F0_ave] = funct_F0_ave(F0)
    Dem = 0 ; 
    F0_ave = 0 ;
for i = 1:length(F0)
    if F0(i) > 0                % Kiem tra F0
        F0_ave = F0_ave + F0(i); % Tinh trung binh F0
        Dem = Dem + 1;            % Tang bien dem de lay gia tri chia
    end
end
    F0_ave = F0_ave/Dem;
end 

