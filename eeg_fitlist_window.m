
% % Load EEG signal

% addpath(genpath('/Users/shadi/Documents/eeg_challenge/training_chall'));
% datalist=dir('/Users/shadi/Documents/eeg_challenge/training_chall');
% no_sub=1;
% datalist = {datalist.name};
% datalist=datalist(1,[7:end]);
% a=[char(datalist(1,no_sub)),'.mat'];
% signals=load(a);
% allsignals=signals.val;
% EEG=allsignals(1:7,:);
% %% arousal lables
% EEG_win=EEG(:,[1:6000]);


function all_feat_mpc_freq=eeg_fitlist_window(EEG_win,fs);

% shape of EEG_win ( no_chanel,segmented EEG) 
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

% the input is a matrix of no_channel*windowsize


no_chan=size(EEG_win,1);
for i=1:no_chan
    for j=1:no_chan
        data1=EEG_win(i,:);
        data2=EEG_win(j,:);
        
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
mpc_delta2=connect_choose(mpc_delta);
mpc_theta2=connect_choose(mpc_theta);
mpc_alpha2=connect_choose(mpc_alpha);
mpc_beta2=connect_choose(mpc_beta);
mpc_gamma2=connect_choose(mpc_gamma);
mpc_all=[mpc_delta2,mpc_theta2,mpc_alpha2,mpc_beta2,mpc_gamma2]; % 21*5=105 features 


    function inp2=connect_choose(inp)
        % order is elements of the lower tringular ( row by row)
        mpc=triu(inp);
        ind=find(mpc(:));
        mpc2=mpc(ind);
        ind2=find(mpc2~=1);
        inp2=mpc2(ind2);
    end






for i=1:no_chan
    %            winval=eeg_filt(i,[startidx:finidx]);
    inp_freq=EEG_win(i,:);
    
    %% EEG power spectrum features
    [s,ff]=pwelch(inp_freq,4*fs,3*fs,[],fs);
    
    pow_delta = bandpower(s,ff,[0.5,4],'psd');
    pow_theta= bandpower(s,ff,[4,8],'psd');
    pow_alpha= bandpower(s,ff,[8,13],'psd');
    pow_beta= bandpower(s,ff,[14,30],'psd');
    pow_gamma= bandpower(s,ff,[32,45],'psd');
    %%freq with respect to max peak
    [max_pow,ind_maxp]=max(s);
    freq_maxpow=ff(ind_maxp);
   fet_eeg_freq{i,1}=[pow_delta,pow_theta,pow_alpha,pow_beta,pow_gamma,freq_maxpow];
    
end
fet_eeg_freq_mat=cell2mat(fet_eeg_freq); % no_chan*6(freq features)
fet_eeg_freq_mat_vec=fet_eeg_freq_mat(:); 

    
    
 
all_feat_mpc_freq=[fet_eeg_freq_mat_vec;mpc_all];
end











