clear all
clc

% cont
load('cont_channels_21.mat')
fs_d=500;
fs=fs_d;
window=blackman(fs_d*4);
overlap=fs_d*3;
mpc_cont_freq=cell(length(cont_channels_21),1);


for k=1:length(cont_channels_21) 
    k
    data=cont_channels_21(k,1);
    data=cell2mat(data);
   
    mpc_all=cell(1,5);
    
for i=1:5%:19 
    for j=1:5%19 
        i
        j
        
        data1=data(i,:);
        data2=data(j,:);
%         
%       [spec_data1,w1,time1]=spectrogram(data1,window,overlap,[],fs_d);
%       [spec_data2,w2,time2]=spectrogram(data2,window,overlap,[],fs_d);
      
      %delta
%        delta1=sig_connect(0.5,4,spec_data1,w1);
%        delta2=sig_connect(0.5,4,spec_data2,w2); 
       
         delta1= eegfilt(data1,fs_d,.5,4);
         delta2= eegfilt(data2,fs,.5,4);


        mpc_delta(i,j) = my_mean_phase_coherence(delta1,delta2);
        
     %theta 
%        theta1=sig_connect(4,8,spec_data1,w1);
%        theta2=sig_connect(4,8,spec_data2,w2); 

        theta1= eegfilt(data1,fs,4,8);
        theta2= eegfilt(data2,fs,.5,4);

        mpc_theta(i,j) = my_mean_phase_coherence(theta1,theta2);
        
         %alpha 
%        alpha1=sig_connect(8,13,spec_data1,w1);
%        alpha2=sig_connect(8,13,spec_data2,w2); 
        alpha1= eegfilt(data1,fs,8,13);
        alpha2= eegfilt(data2,fs,.5,4);

        mpc_alpha(i,j) = my_mean_phase_coherence(alpha1,alpha2);
        
             %beta 
%        beta1=sig_connect(8,13,spec_data1,w1);
%        beta2=sig_connect(8,13,spec_data2,w2); 

        beta1= eegfilt(data1,fs,14,30);
        beta2= eegfilt(data2,fs,14,30);
        
        mpc_beta(i,j) = my_mean_phase_coherence(beta1,beta2);
        
             %gamma
%       gamma1=sig_connect(8,13,spec_data1,w1);
%        gamma2=sig_connect(8,13,spec_data2,w2); 

        gamma1= eegfilt(data1,fs,32,45);
        gamma2= eegfilt(data2,fs,32,45);

        mpc_gamma(i,j)=my_mean_phase_coherence(gamma1,gamma2);
        
        
    end
    mpc_all{1,1}=mpc_delta;
    mpc_all{1,2}=mpc_theta;
    mpc_all{1,3}=mpc_alpha;
   mpc_all{1,4}=mpc_beta;
   mpc_all{1,5}=mpc_gamma;

   



   
end

mpc_cont_freq{k,1}=mpc_all;
end 


% save mpc_cont mpc_cont 