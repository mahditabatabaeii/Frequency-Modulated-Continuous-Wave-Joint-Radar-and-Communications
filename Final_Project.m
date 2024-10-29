clc; clear; close all;

Fs = 2 * 10^9;

IF_bits = 6; 
BW_bits = 7;

f_initial = [(-100)*10^6 : 200*10^6/IF_bits : 100*10^6 - 200*10^6/IF_bits];
Bandwidth = [400*10^6 : 400*10^6/BW_bits : 800*10^6 - 400*10^6/BW_bits];

Tc = 15.2 * 10^(-6); 
T_guard = 0.1 * 10^(-6); 
t = 0 : 1/Fs : 256*Tc + 256*T_guard - 1/Fs;

guard = zeros(1 , T_guard*Fs);
y = [];

Transmitted_chirps=[];
for n = 1:256
    Random_IF = randi(6);
    Random_BW = randi(7);
    Transmitted_chirps=[Transmitted_chirps;f_initial(Random_IF),Bandwidth(Random_BW)];
    f1 = f_initial(1 , Random_IF);
    f2 = f1 + Bandwidth(1 , Random_BW);
    slope = (f2 - f1) / Tc;
    f = f1 : slope/Fs : f2-slope/Fs;
    t_n = t(1 , 1:length(f));
    y_n = 1.*exp(1i.*(2.*pi.*(f1.*t_n+0.5.*((f2-f1)./Tc).*t_n.^2)));
    if n==1
        [ss,ff,tt] = stft(y_n,Fs);
        nexttile
        mesh(tt,ff,abs(ss).^2)
        title("STFT of transmitted chirp with IF="+string(f_initial(Random_IF))+" and BW="+string(Bandwidth(Random_BW)))
        view(2), axis tight
    end
    y = [y y_n guard];
end    


y0 = y + wgn(1,length(y),0,'complex');
y10 = y + wgn(1,length(y),10,'complex');
y20 = y + wgn(1,length(y),20,'complex');

Received_Chirps1=[];
Received_Chirps2=[];
Received_Chirps3=[];
for n = 1:256
    chirp1 = y0((n-1)*(Tc+T_guard)*Fs+1:(n-1)*(Tc+T_guard)*Fs+1+30400-1);
    chirp2 = y10((n-1)*(Tc+T_guard)*Fs+1:(n-1)*(Tc+T_guard)*Fs+1+30400-1);
    chirp3 = y20((n-1)*(Tc+T_guard)*Fs+1:(n-1)*(Tc+T_guard)*Fs+1+30400-1);
    I1=0;
    I2=0;
    I3=0;
    J1=0;
    J2=0;
    J3=0;
    for k = 1:IF_bits
        for j = 1:BW_bits
            f1 = f_initial(1 , k);
            f2 = f1 + Bandwidth(1 , j);
            slope = (f2 - f1) / Tc;
            f = f1 : slope/Fs : f2-slope/Fs;
            t_n = t(1 , 1:length(f));
            y_n = 1.*exp(1i.*(2.*pi.*(f1.*t_n+0.5.*((f2-f1)./Tc).*t_n.^2)));
            if k==1 && j==1
                MSE_chirp1=calculate_mse(chirp1,y_n);
                MSE_chirp2=calculate_mse(chirp2,y_n);
                MSE_chirp3=calculate_mse(chirp3,y_n);
                disp(mse);
                I1=k;
                I2=k;
                I3=k;
                J1=j;
                J2=j;
                J3=j;
            else
                if calculate_mse(chirp1,y_n)<MSE_chirp1
                    MSE_chirp1=calculate_mse(chirp1,y_n);
                    I1=k;
                    J1=j;
                end
                if calculate_mse(chirp2,y_n)<MSE_chirp2
                    MSE_chirp2=calculate_mse(chirp2,y_n);
                    I2=k;
                    J2=j;
                end
                if calculate_mse(chirp3,y_n)<MSE_chirp3
                    MSE_chirp3=calculate_mse(chirp3,y_n);
                    I3=k;
                    J3=j;
                end
            end
        end
    end
    if n==1
        [ss1,ff1,tt1] = stft(chirp1,Fs);
        nexttile
        mesh(tt1,ff1,abs(ss1).^2)
        title("STFT of received chirp with 0dB SNR")
        view(2), axis tight
        
        [ss2,ff2,tt2] = stft(chirp2,Fs);
        nexttile
        mesh(tt2,ff2,abs(ss2).^2)
        title("STFT of received chirp with 10dB SNR")
        view(2), axis tight
        
        [ss3,ff3,tt3] = stft(chirp3,Fs);
        nexttile
        mesh(tt3,ff3,abs(ss3).^2)
        title("STFT of received chirp with 20dB SNR")
        view(2), axis tight
    end
    Received_Chirps1=[Received_Chirps1;f_initial(I1),Bandwidth(J1)];
    Received_Chirps2=[Received_Chirps2;f_initial(I2),Bandwidth(J2)];
    Received_Chirps3=[Received_Chirps3;f_initial(I3),Bandwidth(J3)];
end

function mse = calculate_mse(x, y)
diff = x - y;
squared_diff = diff .* conj(diff);
mse = mean(squared_diff);
end

    