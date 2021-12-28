%This program simulates the working of the ASK-modulated digital communication system
%****************************************************************************
%The maximum bandwidth of the message signal is assumed to be fm=50 Hz
%The message signal is sampled at a sampling rate of fs=10*fm
%The samples are encoded in binary and Hamming parity bits are generated
%The binary Hamming codes are encoded using NRZ on-off line coding scheme
%The encoded samples are then modulated using ASK
%Additive White Gaussian Noise is assumed to be present in the channel with a given SNR
%The reciever rectifies the signal and uses envelope detection to detect the on-off signal
%The encoded signal is checked for errors using Hamming parity bits
%The error-corrected encoded samples are decoded to decimal values
%The sampled signal is generated from the above decimal values
%The sampled signal is passed through a lowpass filter to obtain the reconstructed signal
%Octave packages signal and communications have been used
fm=50;
Tm=1/fm;
fs=10*fm;
Ts=1/fs;
t=0:Ts/70:2*Tm; %creating the time vector
m=6+5*cos(2*pi*30*(t-0.02))+sin(2*pi*20*t); %message signal to be transmitted
fc=10000;
l=length(m);
sampled_m=zeros(1,l);
samples=zeros(1,((2*Tm)/Ts)+1);
i=2;
sampled_m(1)=m(1);
while(i<=l)
if(rem(i-1,70)==0)
sampled_m(i)=m(i);
samples((i-1)/70)=m(i);
else
sampled_m(i)=sampled_m(i-1);
endif;
i=i+1;
endwhile; %sampled signal generated
sampled_m_bin=dec2bin(floor(samples)); %binary codes for the samples obtained
p1=xor(str2num(sampled_m_bin(:,4)),xor(str2num(sampled_m_bin(:,3)),str2num(sampled_m_bin(:,1))));
p2=xor(str2num(sampled_m_bin(:,4)),xor(str2num(sampled_m_bin(:,2)),str2num(sampled_m_bin(:,1))));
p4=xor(str2num(sampled_m_bin(:,3)),xor(str2num(sampled_m_bin(:,2)),str2num(sampled_m_bin(:,1))));
d3=str2num(sampled_m_bin(:,4));
d5=str2num(sampled_m_bin(:,3));
d6=str2num(sampled_m_bin(:,2));
d7=str2num(sampled_m_bin(:,1));
encoded_samples=[d7,d6,d5,p4,d3,p2,p1]; %Hamming codes generated for the samples
pcm=zeros(1,l);
pcm(1)=encoded_samples(1,1);
encoded_sample=encoded_samples(1,:);
j=2;k1=1;
while(j<=l)
if(rem(j-1,70)==0)
encoded_sample=encoded_samples((j-1)/70,:); k1=1;
endif
if(rem(j-1,10)==0)
pcm(j)=encoded_sample(k1);
k1=k1+1;
else
pcm(j)=pcm(j-1);
endif;
j=j+1;
endwhile; %NRZ on-off line codes obtained
ask_output=sin(2*pi*fc*t).*(pcm>0); %ASK modulated signal transmitted
%end of modulation
snr=50;
recieved=awgn(ask_output,snr); %AWGN corrupts signal
th=0.2; %threshold level for detection
rect_recieved=recieved.*(recieved>0); %recieved signal rectified
env=abs(hilbert(rect_recieved)); %envelope detection used
rencoded_samples=zeros(((2*Tm)/Ts+1),7);
k3=1;
k=2;
i1=1;
while(k<=l)
if(rem((k-1),5)==0 && rem((k-1),10)!=0)
rencoded_samples(k3,i1)=(env(k)>th);
i1=i1+1;
endif;
if(rem((k-1),70)==0)
i1=1;k3=k3+1;
endif;
k=k+1;
endwhile;
rp4=xor(rencoded_samples(:,1),rencoded_samples(:,2),rencoded_samples(:,3),rencoded_samples(:,4));
rp2=xor(rencoded_samples(:,1),rencoded_samples(:,2),rencoded_samples(:,5),rencoded_samples(:,6));
rp1=xor(rencoded_samples(:,1),rencoded_samples(:,3),rencoded_samples(:,5),rencoded_samples(:,7));
errpos=bin2dec(num2str([rp4,rp2,rp1]));
c=(2*Tm)/Ts+1;
j1=1;
while(j1<=c)
if(errpos(j1)==3||errpos(j1)==5||errpos(j1)==6||errpos(j1)==7)
rencoded_samples(j1,8-errpos(j1))=not(rencoded_samples(j1,8-errpos(j1)));
endif;
j1=j1+1;
endwhile; %error-corrected encoded samples
rsamples=bin2dec(num2str([rencoded_samples(:,1),rencoded_samples(:,2),rencoded_samples(:,3),rencoded_samples(:,5)]));
i=2;
rsampled=zeros(1,l);
rsampled(1)=rsamples(1);
while(i<=l)
if(rem(i-1,70)==0)
rsampled(i)=rsamples((i-1)/70);
else
rsampled(i)=rsampled(i-1);
endif;
i=i+1;
endwhile; %sampled signal at the reciever obtained from error-corrected encoded samples
Lfft=2^ceil(log2(l)+1);
H=zeros(1,Lfft);
H(Lfft/2-fm:Lfft/2+fm)=1;
SGH=fftshift(fft(rsampled,Lfft));
reconstructed_M=SGH.*H;
reconstructed_m=real(ifft(fftshift(reconstructed_M)));
reconstructed_m=reconstructed_m(1:l);
reconstructed_m=abs(hilbert(reconstructed_m)); %signal reconstructed using lowpass filter
%end of demodulation
figure(1);
subplot(221); td1=plot(t,m);
xlabel('\it t(in s)');
ylabel('m(t)');
title('Original message signal');
subplot(222); td2=plot(t,sampled_m);
xlabel('\it t(in s)');
ylabel('sampled_m(t)');
title('Sampled message signal');
subplot(224); td3=plot(t,rsampled);
xlabel('\it t(in s)');
ylabel('rsampled(t)');
title('Sampled signal at the reciever');
subplot(223); td4=plot(t,reconstructed_m);
xlabel('\it t(in s)');
ylabel('reconstructed_m(t)');
title('Reconstructed signal at the reciever');
figure(2)
plot(t,pcm);
xlabel('\it t(in s)');
ylabel('pcm(t)');
title('PCM NRZ on-off signal generated');
figure(3)
plot(t,ask_output);
xlabel('\it t(in s)');
ylabel('ask_output(t)');
title('ASK modulated transmitted signal');
figure(4)
plot(t,recieved);
xlabel('\it t(in s)');
ylabel('recieved(t)');
title('Recieved signal corrupted by AWGN');
figure(5)
plot(t,env);
xlabel('\it t(in s)');
ylabel('env(t)');
title('Envelope detected signal');
