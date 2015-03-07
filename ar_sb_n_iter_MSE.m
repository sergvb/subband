% get access to model
%clc, clear all, clf ;
curPath = pwd() ;
cd('../phd_src-master/src/tests/tsim/model') ;
modelPath = pwd() ;
cd( curPath ) ;
addpath(modelPath) ;

num_of_tests = 1000 ;
N = 16 ;
SNR_dB = -17:1:14 ;

fd = 16E3 ;
fs = 2.3E3 ;
w1 = -2*pi*4E3/fd;
w2 = 2*pi*4E3/fd;

A = 1 ; E = A^2 / 2 ;
sigma = E ./ (10 .^ (SNR_dB./10)) ;
SNR = E ./ sigma ;

% phase_arg = 2*pi*1*fs/fd*(0:N-1) ;
% s = A * cos(phase_arg) ;

freq1 = zeros(length(SNR_dB), 1) ;
freq2 = zeros(length(SNR_dB), 1) ;
freq3 = zeros(length(SNR_dB), 1) ;

freq1_4N = zeros(length(SNR_dB), 1) ;
freq2_4N = zeros(length(SNR_dB), 1) ;
freq3_4N = zeros(length(SNR_dB), 1) ;

freq4 = zeros(length(SNR_dB), 1) ;



EE = zeros(N, N, N);
for tau = 0:N-1
    for k = 1:N
        for n = 1:N
            d = k - n + tau;
            EE(k, n, tau + 1) = 1.3;
            if d ~= 0
                EE(k, n, tau + 1) = -1i/d*(exp(1i*w2*d) - exp(1i*w1*d));
            end
        end
    end
end

pObj = parpool(7) ;

parfor jj=1:length(SNR_dB)
    E = EE;
    r = zeros(1, N);
    
    fprintf('Actual: %.4f  dB\n', SNR_dB(jj));
    
    for k=1:num_of_tests
        fs = 2000 + 1000*rand(1) ;
        phase_arg = 2*pi*1*fs/fd*(0:N-1) ;
        s = A * cos(phase_arg) ;
        x = s + sqrt(sigma(jj))*(randn(size(s))) ;
        
        %%%%%%%%%%%%%%%%%%%
        % signal + noise
        X = fft(x) ;
        XX = X.*conj(X) ;
        rxx1 = ifft(XX) ;
        rxx2 = ifft(XX .^ 4) ;
        rxx3 = ifft(XX .^ 8) ;
                
        %%%%%%%%%%%%
        % 1
        b1 = ar_model([rxx1(1); rxx1(2); rxx1(3)]) ;
        [poles1, omega0_1, Hjw0_1] = get_ar_pole(b1) ;
        freq1(jj) = freq1(jj) + (omega0_1*fd/2/pi - fs)^2 ;
        
        %%%%%%%%%%%%
        % 2
        b2 = ar_model([rxx2(1); rxx2(2); rxx2(3)]) ;
        [poles2, omega0_2, Hjw0_2] = get_ar_pole(b2) ;
        freq2(jj) = freq2(jj) + (omega0_2*fd/2/pi - fs)^2 ;
        
        %%%%%%%%%%%%
        % 3
        b3 = ar_model([rxx3(1); rxx3(2); rxx3(3)]) ;
        [poles3, omega0_3, Hjw0_3] = get_ar_pole(b3) ;
        freq3(jj) = freq3(jj) + (omega0_3*fd/2/pi - fs)^2 ;
        
        
        %%%%%%%%%%%%
        % 4*N approx
        X_4N = fft(x, 4*N) ;
        XX_4N = X_4N.*conj(X_4N) ;
        rxx1_4N = ifft(XX_4N) ;
        rxx2_4N = ifft(XX_4N .^ 4) ;
        rxx3_4N = ifft(XX_4N .^ 8) ;
        
        %%%%%%%%%%%%
        % 1_4N
        b1_4N = ar_model([rxx1_4N(1); rxx1_4N(2); rxx1_4N(3)]) ;
        [poles1_4N, omega0_1_4N, Hjw0_1_4N] = get_ar_pole(b1_4N) ;
        freq1_4N(jj) = freq1_4N(jj) + (omega0_1_4N*fd/2/pi - fs)^2 ;
        
        %%%%%%%%%%%%
        % 2_4N
        b2_4N = ar_model([rxx2_4N(1); rxx2_4N(2); rxx2_4N(3)]) ;
        [poles2_4N, omega0_2_4N, Hjw0_2_4N] = get_ar_pole(b2_4N) ;
        freq2_4N(jj) = freq2_4N(jj) + (omega0_2_4N*fd/2/pi - fs)^2 ;
        
        %%%%%%%%%%%%
        % 3_4N
        b3_4N = ar_model([rxx3_4N(1); rxx3_4N(2); rxx3_4N(3)]) ;
        [poles3_4N, omega0_3_4N, Hjw0_3_4N] = get_ar_pole(b3_4N) ;
        freq3_4N(jj) = freq3_4N(jj) + (omega0_3_4N*fd/2/pi - fs)^2 ;
        
        
        %%%%%%%%%%%%
        % subband
        for ccc = 1:7
            for tau = 1:N
                r(tau) = 1/(2*pi)*x*E(:,:,tau)*x.';
            end
            x = r;
        end
        
        rxx4 = [ r(1); r(2); r(3); ];
        
        b4 = ar_model([rxx4(1); rxx4(2); rxx4(3)]) ;
        [poles4, omega0_4, Hjw0_4] = get_ar_pole(b4) ;
        freq4(jj) = freq4(jj) + (omega0_4*fd/2/pi - fs)^2 ;
        
        %fprintf('Estimated freq: %.4f Hz\n', freq3(k));
    end ;
    
    freq1(jj) = sqrt(freq1(jj) / num_of_tests) ;
    freq2(jj) = sqrt(freq2(jj) / num_of_tests) ;
    freq3(jj) = sqrt(freq3(jj) / num_of_tests) ;
    
    freq1_4N(jj) = sqrt(freq1_4N(jj) / num_of_tests) ;
    freq2_4N(jj) = sqrt(freq2_4N(jj) / num_of_tests) ;
    freq3_4N(jj) = sqrt(freq3_4N(jj) / num_of_tests) ;
    
    freq4(jj) = sqrt(freq4(jj) / num_of_tests) ;
    
end ; % SNR

delete(pObj);

SNR_dB_acf = SNR_dB ;
save('freq_sko_ar', 'freq1', 'freq2', 'freq3', ...
    'freq1_4N', 'freq2_4N', 'freq3_4N', ...
    'freq4', 'SNR_dB_acf')

figure(1) ,
semilogy(SNR_dB, freq1, '-go', ...
    SNR_dB, freq2, '-b*', ...
    SNR_dB, freq3, '-r+', ...
    SNR_dB, freq1_4N, '-c.', ...
    SNR_dB, freq2_4N, '-kd', ...
    SNR_dB, freq3_4N, '-y<', ...
    SNR_dB, freq4, '-m^') ,
title('MSE') ,
legend('1', '2', '3', '1_{4N}', '2_{4N}', '3_{4N}', 'SB') ;
grid on;

% remove model path
rmpath(modelPath) ;