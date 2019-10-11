clear all
clc

addpath(genpath('/Users/shadi/Documents/eeg_challenge/training_chall'));

datalist=dir('/Users/shadi/Documents/eeg_challenge/training_chall');
datalist = {datalist.name};
datalist=datalist(1,[7:end]);
%% arousal lables


fs=200;
windowSize=30;
win_size_fs=windowSize*fs;

no_chan=7;


features_eeg_win_all=cell(size(datalist,2),1);
features_eeg_all_chan_connect=cell(size(datalist,2),1);
for no_sub=1:size(datalist,2)
    no_sub
    a=[char(datalist(1,no_sub)),'.mat'];
    signals=load(a);
    allsignals=signals.val;
    EEG=allsignals(1:7,:);
    lenght=size(EEG,2);
    
    
    %% preprocessing for the whole time winfow
    
    % filtering
    [filt] = lowpassfilter(EEG,fs,45);
    [eeg_filt] = highpassfilter(filt,fs,0.5);
    %     % rejecting low quality ? ask Albe ( from Nadi's paper)
    %     ff=diff(eeg_filt);
    %     test=eeg_filt(1,:);
    %
    %     [ind,vv]=find(ff(1,:)>100);
    %     test(vv)=0;
    %     figure
    %     plot(eeg_filt(1,:))
    %     hold on
    %     plot(test(1,:),'r')
    
    % other filtering codes
    % [smoothdata,filtwts] = eegfilt_v(EEG,fs,0.5,45);
    
    % cutoff_low=45;
    % cutoff_high=0.5;
    % order=2;
    % [B_low,A_low] = butter(order,2*cutoff_low/fs,'low');
    % low_passed = filtfilt(B_low,A_low,EEG);
    %
    % [B_high,A_high] = butter(order,2*cutoff_high/fs,'high');
    % filt_data = filtfilt(B_high,A_high,low_passed);
    
    
 
    k=1;
    features_eeg_win=[];
    features_eeg_chan_connect=[];
    features_eeg=cell(no_chan,1);
    
    startidx=1;
    finidx=win_size_fs;
    deltawindow=win_size_fs/2; % in case of overlap
    
    
    
    
    
    
    
    
target=[];
    while finidx<=lenght
        k
        % channel orders: 'F3-M2'	'F4-M1'	'C3-M2'	'C4-M1'	'O1-M2'	'O2-M1'	'E1-M2'
        
         %%finding target arousals in specified windows
            arous=[char(datalist(1,no_sub)),'-arousal.mat'];
            arous=load(arous);
            arous_label=arous.data.arousals;
            currLabelVec=arous_label(1,[startidx:finidx]);
            arous=find(currLabelVec==1);
            undefined=find(currLabelVec==-1);
            non_arous=find(currLabelVec==0);
            lable_vote=[sum(currLabelVec==-1),sum(currLabelVec==0),sum(currLabelVec==1)];
            idx=argmax(lable_vote);
            idx=idx-1;  % idx=0 : undefined / idx=1 non_aours/ idx=1 arous
            target=[target,idx];
        
        
        %% connectivity for 5 bands ( 7*7 matrix) (7:no_channel)
        winval_all=eeg_filt(:,[startidx:finidx]);  %all 7 channels together 
        for i=1:no_chan
            for j=1:no_chan
                data1=winval_all(i,:);
                data2=winval_all(j,:);
                
                % delta
                
                [data1_low] = lowpassfilter(data1,fs,4);
                [delta1] = highpassfilter(data1_low,fs,0.5);
                
                [data2_low] = lowpassfilter(data2,fs,4);
                [delta2] = highpassfilter(data2_low,fs,0.5);
                % theta
                
                [data1_low] = lowpassfilter(data1,fs,8);
                [theta1] = highpassfilter(data1_low,fs,4);
                
                [data2_low] = lowpassfilter(data2,fs,8);
                [theta2] = highpassfilter(data2_low,fs,4);
                % alpha
                
                [data1_low] = lowpassfilter(data1,fs,13);
                [alpha1] = highpassfilter(data1_low,fs,8);
                
                [data2_low] = lowpassfilter(data2,fs,13);
                [alpha2] = highpassfilter(data2_low,fs,8);
                % beta
                
                [data1_low] = lowpassfilter(data1,fs,30);
                [beta1] = highpassfilter(data1_low,fs,14);
                
                [data2_low] = lowpassfilter(data2,fs,30);
                [beta2] = highpassfilter(data2_low,fs,14);
                % gamma
                
                [data1_low] = lowpassfilter(data1,fs,45);
                [gamma1] = highpassfilter(data1_low,fs,32);
                
                [data2_low] = lowpassfilter(data2,fs,45);
                [gamma2] = highpassfilter(data2_low,fs,32);
                
                
                
                
                mpc_delta(i,j) = my_mean_phase_coherence(delta1,delta2);
                mpc_theta(i,j) = my_mean_phase_coherence(theta1,theta2);
                mpc_alpha(i,j) = my_mean_phase_coherence(alpha1,alpha2);
                mpc_beta(i,j) = my_mean_phase_coherence(beta1,beta2);
                mpc_gamma(i,j) = my_mean_phase_coherence(gamma1,gamma2);
                
                
                
                
            end
        end
        features_eeg_chan_connect{k,1}=(mpc_delta);
        features_eeg_chan_connect{k,2}=(mpc_theta);
        features_eeg_chan_connect{k,3}=(mpc_alpha);
        features_eeg_chan_connect{k,4}=(mpc_beta);
        features_eeg_chan_connect{k,5}=(mpc_gamma);
        
        
        
        
        
        for i=1:no_chan
            winval=eeg_filt(i,[startidx:finidx]);
            
           
            %% EEG power spectrum features
            [s,ff]=pwelch(winval,4*fs,3*fs,[],fs);
            
            pow_delta = bandpower(s,ff,[0.5,4],'psd');
            pow_theta= bandpower(s,ff,[4,8],'psd');
            pow_alpha= bandpower(s,ff,[8,13],'psd');
            pow_beta= bandpower(s,ff,[14,30],'psd');
            pow_gamma= bandpower(s,ff,[32,45],'psd');
            %%freq with respect to max peak
            [max_pow,ind_maxp]=max(s);
            freq_maxpow=ff(ind_maxp);
            features_eeg{i,1}=[pow_delta,pow_theta,pow_alpha,pow_beta,pow_gamma,freq_maxpow];
            features_eeg_mat=cell2mat(features_eeg);
            
            
        end
        features_eeg_win{k,1}=features_eeg_mat;
        
        startidx=startidx+deltawindow;
        finidx=finidx+deltawindow;
        k=k+1;
        
    end
    
    features_eeg_win_all{no_sub,1}=features_eeg_win;
    features_eeg_win_all{no_sub,2}=features_eeg_chan_connect;
    target_all{no_sub,1}=target;
end


%% rearranging for classifier 
%  output: a cell of nosub*[(no_win)*(147 features)]
load('features_eeg_win_all.mat')
featurelist={'PowDelta_c1','PowTetha_c1','PowAlpha_c1','PowBeta_c1','PowGamma_c1','freq_c1',...
    'PowDelta_c2','PowTetha_c2','PowAlpha_c2','PowBeta_c2','PowGamma_c2','freq_c2',...
    'PowDelta_c3','PowTetha_c3','PowAlpha_c3','PowBeta_c3','PowGamma_c3','freq_c3',...
    'PowDelta_c4','PowTetha_c4','PowAlpha_c4','PowBeta_c4','PowGamma_c4','freq_c4',...
    'PowDelta_c5','PowTetha_c5','PowAlpha_c5','PowBeta_c5','PowGamma_c5','freq_c5',...
    'PowDelta_c6','PowTetha_c6','PowAlpha_c6','PowBeta_c6','PowGamma_c6','freq_c6',...
    'PowDelta_c7','PowTetha_c7','PowAlpha_c7','PowBeta_c7','PowGamma_c7','freq_c7',...
    'connect(21)(delta)','connect(21)(delta)','connect(21)(delta)','connect(21)(delta)', 'connect(21)(delta)','connect(21)(delta)','connect(21)(delta)','connect(21)(delta)', 'connect(21)(delta)','connect(21)(delta)','connect(21)(delta)','connect(21)(delta)', 'connect(21)(delta)','connect(21)(delta)','connect(21)(delta)','connect(21)(delta)', 'connect(21)(delta)','connect(21)(delta)','connect(21)(delta)','connect(21)(delta)','connect(21)(delta)',...
    'connect(21)(theta)','connect(21)(theta)','connect(21)(theta)','connect(21)(theta)', 'connect(21)(theta)','connect(21)(theta)','connect(21)(theta)','connect(21)(theta)', 'connect(21)(theta)','connect(21)(theta)','connect(21)(theta)','connect(21)(theta)', 'connect(21)(theta)','connect(21)(theta)','connect(21)(theta)','connect(21)(theta)', 'connect(21)(theta)','connect(21)(theta)','connect(21)(theta)','connect(21)(theta)','connect(21)(theta)',...
    'connect(21)(alpha)','connect(21)(alpha)','connect(21)(alpha)','connect(21)(alpha)', 'connect(21)(alpha)','connect(21)(alpha)','connect(21)(alpha)','connect(21)(alpha)', 'connect(21)(alpha)','connect(21)(alpha)','connect(21)(alpha)','connect(21)(alpha)', 'connect(21)(alpha)','connect(21)(alpha)','connect(21)(alpha)','connect(21)(alpha)', 'connect(21)(alpha)','connect(21)(alpha)','connect(21)(alpha)','connect(21)(alpha)','connect(21)(alpha)',...
    'connect(21)(beta)','connect(21)(beta)','connect(21)(beta)','connect(21)(beta)', 'connect(21)(beta)','connect(21)(beta)','connect(21)(beta)','connect(21)(beta)', 'connect(21)(beta)','connect(21)(beta)','connect(21)(beta)','connect(21)(beta)', 'connect(21)(beta)','connect(21)(beta)','connect(21)(beta)','connect(21)(beta)', 'connect(21)(beta)','connect(21)(beta)','connect(21)(beta)','connect(21)(beta)','connect(21)(beta)',...  
'connect(21)(gamma)','connect(21)(gamma)','connect(21)(gamma)','connect(21)(gamma)', 'connect(21)(gamma)','connect(21)(gamma)','connect(21)(gamma)','connect(21)(gamma)', 'connect(21)(gamma)','connect(21)(gamma)','connect(21)(gamma)','connect(21)(gamma)', 'connect(21)(gamma)','connect(21)(gamma)','connect(21)(gamma)','connect(21)(gamma)', 'connect(21)(gamma)','connect(21)(gamma)','connect(21)(gamma)','connect(21)(gamma)','connect(21)(gamma)'};


% no: 6*7+21*5=147 features 
% (1,2)to(1,7),(2,3)to(2,7),(3,4)to(3,7),(4,5)to(4,7),(5,6)to(5,7),(6,7)



% no_sub=size(features_eeg_win_all,1);
no_sub=2;
matr_all=cell(no_sub,1);

for mm=1:no_sub
    
    matr1=[];
matr2=[];
    
for i=1:length(features_eeg_win_all{mm,1})
    
    
    %% psd features 
    matr1=[matr1;features_eeg_win_all{mm,1}{i}(:)'];
    %% connectivity features 
    sel_fin_all=cell(1,5);
    for j=1:5 % 5 bands for connectivity 
    
    sel=triu(features_eeg_win_all{mm,2}{i,j});
    sel=sel(:);
    sel2=find(sel);
    sel3=sel(sel2);
    sel4=find(sel3~=1);
    sel_fin=sel3(sel4);
    sel_fin_all{1,j}=sel_fin;
    end
    sel_fin_all_mat=cell2mat(sel_fin_all');
    
   matr2=[matr2;sel_fin_all_mat'];
   matr_tog=[matr1,matr2];
   
end
   matr_all{mm,1}=matr_tog;
end

% matr_all :cell(nosub,feature_matrix) , feature_matrix( no_wind, 147
% features),

% save target_all target_all
%         
%     end

    
% EEG_f_all.feat=featurelist';

    
    
    


