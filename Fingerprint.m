%%
%% =====================================================================================
%%       Filename:  Fingerprint.m 
%%
%%    Description:  Micro-CSI extraction
%%
%%         Author:  Ruiqi Kong 
%%         Email :  <kr020@ie.cuhk.edu.hk>
%%   Organization:  WiNS group @ The chiniese university of hong kong
%%
%%   Copyright (c)  WiNS group @ The chiniese university of hong kong
%% =====================================================================================
%%

classdef Fingerprint < handle

    properties
        N_csi = 1;
        devices=struct;
        devices_list=0
        usedrx = 1;
        oe = 1;
        A = [];
    end
    
    methods
        function obj = Fingerprint(N_csi,usedrx,oe,n_taps)
            obj.N_csi = N_csi;
            obj.usedrx = usedrx;
            obj.oe = oe;
            DFT=dftmtx(64);
            scn = 52;
            F=DFT([(64-scn/2+1):64,2:(1+scn/2)],[1:n_taps+1,64-n_taps+1:64]);
            obj.A=F*inv(F'*F)*F';
        end

        function get_micro_csi_group(obj,CSI)
            microcsi_group={};
            for n=1:size(CSI,2)
                microcsi = get_micro_csi(obj,CSI{1,n});
                microcsi_group=cat(2,microcsi_group,microcsi);
            end
            pre_device=obj.devices_list(end);
            obj.devices.("device"+string(pre_device+1))={microcsi_group};
            obj.devices_list=1:pre_device+1;       
        end

        function microcsi = get_micro_csi(obj,CSI)
            %input csi data   [num_sc, num_tx, num_rx, num_csi];        
            %reshape csi format
            est=permute(CSI, [4 1 3 2]);  % [num_csi, num_sc, num_rx, num_tx];

            for ntx=1:size(est,4)
                f=microcsi_extraction(obj,squeeze(est(:,:,:,ntx)),obj.N_csi);
                microcsi=permute(f,[1 3 2]);
            end
            microcsi=reshape(microcsi,[],ntx,size(f,3),size(f,2));
        end

        function f=microcsi_extraction(obj,csi,num)
                % reshape data as N_csi as one group
                csi_groups=reshape(csi(1:floor(size(csi,1)/num)*num,:,:),num,[],52,length(obj.usedrx));
                % N_m = N_csi * N_rx
                csi_groups=reshape(permute(csi_groups,[4 1 2 3]),num*length(obj.usedrx),[],52,1);           
                f=[];
                for group=1:size(csi_groups,2)
                    csi=squeeze(csi_groups(:,group,:));
                    H_ls=(obj.A*csi.').';
                    microcsi=csi./H_ls;
                    if obj.oe==1
                        gphase=gradient(microcsi);
                        variance=var(gphase.');
                        index=find(variance<2*10^(-3));
                        
                        if isempty(index)
                            % no accurate micro-CSI received. 
                            % require retransmission                         
                            continue;
                        else
                            microcsi=microcsi(index,:);
                            z=zscore((microcsi));microcsi(abs(z)>1)=nan;
                            microcsi=mean(microcsi,1,'omitnan');
                        end
                    else
                        microcsi=mean(microcsi,1,'omitnan');
                    end
                    f=cat(1,f,microcsi);
               end
          end
    end
end