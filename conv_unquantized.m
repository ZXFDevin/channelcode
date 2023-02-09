% BPSK modulation and Convolutional Coding [171 133]
% this can serve as a template for 1/2 rate coding

clear all;

Eb=2; %energy per bit, power=1, 1/2 code rate, Eb=2

SNR=[2,3,4];
loop_max=20000; %total frames
N=314;  %data frame length 
NN=320*2;  %coded frame length
TR = poly2trellis(7,[171 133]); % Define trellis

for loop_snr=1:length(SNR)
    snr=10^(SNR(loop_snr)/10);           %snr in ratio.snr=Eb/N0;
    sigma=sqrt(Eb/(2*snr));              % N0=2*sigma^2£¬Eb/N0=snr
    frame_error=0; 

    for loop_frame=1:loop_max

        source=randi(2,N,1)-1;             %data source
        source0=[source;0;0;0;0;0;0];    %zero padding
                                
        code0=convenc(source0,TR);      %encode frame by frame

        R=1-2*code0;           %BPSK 
                                   
        noise=sigma*(randn(NN,1)+i*randn(NN,1));  %complex AWGN    
        RR=R+noise;                                   %add AWGN noise
   
        RRR=real(RR);  %in BPSK , the receiver just process real part, imag part is all noise.
   
        decoded = vitdec(RRR,TR,35,'term','unquant'); % Decode frame by frame
    
        number0 = sum(abs(decoded(1:N)-source));
        if number0>0
            frame_error=frame_error+1;              
        end

  
        if (mod(loop_frame,2000)==0)
            SNR(loop_snr)
            loop_frame
        end
        
    end

    fer(loop_snr)=frame_error/loop_max;
   
end 

% plot results
semilogy(SNR,fer,'b');
grid on