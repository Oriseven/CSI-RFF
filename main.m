%%
%% ===================================================================================== 
%%       Filename:  main.m 
%%
%%    Description:  CSI-RFF implementation for 'iforest','lof','ocsvm','knn' algorithms
%%
%%         Author:  Ruiqi Kong 
%%         Email :  <kr020@ie.cuhk.edu.hk>
%%   Organization:  WiNS group @ The chiniese university of hong kong
%%
%%   Copyright (c)  WiNS group @ The chiniese university of hong kong
%% =====================================================================================
%%
%% dataloader
clear;
load("CSI_data.mat");
% -------------------
% Each row represent CSI from one NIC; Each column represent CSI collected
% in one condition.
% NICs_order =["ESP32C1","ESP32C2","ESP32C3","ESP32C4","ESP32C5",...
%     "AX200C1","AX200C2","AC8260C1","AC7260C1",...
%     "AC7265C1","RTL8812BU","AR9271C1","AR9271C2","AR9271C3","AR9271C4"];
% Conditions_order =
% ["RoomA_static","RoomA_static","RoomA_mobile","RoomA_mobile","RoomB_static","RoomB_static","RoomB_mobile","RoomB_mobile"];
%% fingerprint construction
% Initialize CSI fingerprint extractor
N_csi = 20; % the number of CSI measurements used for fingerprint construction
N_rx = 1:4; % used rx chains
enable_oe = 1; % enable outlier elimination. (Algorithm 1 in CSI-RFF paper)
n_taps = 8; % the number leakaged taps caused by pulse shaping
fingerprints=Fingerprint(N_csi,N_rx,enable_oe,n_taps);
for nic=1:size(CSI,1) 
    get_micro_csi_group(fingerprints,CSI(nic,:));
end
clearvars -except fingerprints;
%% fingerprint normalization
data=struct2cell(fingerprints.devices);
for i=1:length(data)
    for j= 1:length(data{i,1}{1,1})
        data{i,1}{1,1}{1,j}=zscore((data{i,1}{1,1}{1,j}),[],4);
    end
end
clearvars -except fingerprints data;
%% authentication score
distance={'Euclidean_distance','Manhattan_distance','Chebyshev_distance','Euclidean_angle','Hermitian_angle'};

test_enviroment = [5,6;7,8].'; % roomB
scores_knn={};scores_iforest={};scores_lof={};scores_ocsvm={};
for env = 1:size(test_enviroment,2)
    for dis = 2
        for legal = 1: length(data)
            train_xdata=[];train_ylabel=[];
            f=squeeze(cell2mat(data{legal,1}{1,1}(1,[1:4]).'));% use roomA data for training
            train_xdata=cat(1,train_xdata,f);
            train_xdata=cat(2,real(train_xdata),imag(train_xdata));
            for test_device = 1: length(data)
                test_xdata=[];test_ylabel=[];
                f=squeeze(cell2mat(data{test_device,1}{1,1}(1,test_enviroment(:,env)).'));
                test_xdata=cat(1,test_xdata,f);
                test_xdata=cat(2,real(test_xdata),imag(test_xdata));

            if dis==1
               scores_iforest{dis,env,test_device,legal}=novelty_detection(train_xdata,test_xdata,'iforest');
            end
            scores_lof{dis,env,test_device,legal}=novelty_detection(train_xdata,test_xdata,'lof',dis);
            
            scores_ocsvm{dis,env,test_device,legal}=novelty_detection(train_xdata,test_xdata,'ocsvm',dis);
             
            scores_knn{dis,env,test_device,legal}=novelty_detection(train_xdata,test_xdata,'knn',dis);
            end
        end
    end
end
clearvars -except fingerprints data scores_iforest scores_lof scores_ocsvm scores_knn;
%% authentication results
distance={'Euclidean_distance','Manhattan_distance','Chebyshev_distance','Euclidean_angle','Hermitian_angle'};
ND = 'knn'; %'iforest','lof','ocsvm','knn'
op_far=1; % set false alarm rate / false reject rate; op_far = 0 is not support, please replace 0 with a small number like 0.0001
for env = 1:2 % 1 for static; 2 for  mobile
    disp(['-------------env:   ', num2str(env), '--------------']);
    for dis = 2
        eval(['scores = scores_' ND ';']);
        EER=[];adr=[];nc=1:15;
        for legal_nic=nc
            legitimate=scores{dis,env,legal_nic,legal_nic};
            for i_nic=setdiff(nc,legal_nic)   
            attack=scores{dis,env,i_nic,legal_nic};
            adr(end+1) = adr_calculate(legitimate,attack,op_far,10000);
            end
        end
        ADR_summary(dis,env,:)=[mean(adr(:)),max(adr(:)),min(adr(:))];
        disp(['---------alg: ' ND '----distance: ', distance{dis}, '---env:   ', num2str(env), '----------']);
        disp(['When FRR<=' num2str(op_far) '%, average ADR= ' num2str(ADR_summary(dis,env,1)) '  Max ADR= ' num2str(ADR_summary(dis,env,2)) ' MinADR= ' num2str(ADR_summary(dis,env,3))])
    end
end
