%%
%% =====================================================================================
%%       Filename:  adr_calculate.m 
%%
%%    Description:  attack detection rate calculation
%%
%%         Author:  Ruiqi Kong 
%%         Email :  <kr020@ie.cuhk.edu.hk>
%%   Organization:  WiNS group @ The chiniese university of hong kong
%%
%%   Copyright (c)  WiNS group @ The chiniese university of hong kong
%% =====================================================================================
%%

function adr=adr_calculate(legitimates,attackers,OPvalue,pas0)

m0 = min(legitimates)-1;
num_legitimates = length(legitimates);
m1 = max(attackers)+1;
num_attackers = length (attackers);

% calculation of the step
pas1 = (m1 - m0)/pas0;
x = [m0:pas1:m1].';
num = length (x);

%-------- calculation of FAR and FRR
FRR=100;FAR=100;
if num~=0
    for i=1:num
        fr=0;
        fa=0;
        for j=1:num_legitimates
            if legitimates(j)>x(i)
                fr=fr+1;
            end
        end
        for k=1:num_attackers
            if attackers(k)<=x(i)
                fa=fa+1;
            end
        end
        FRR(i)=100*fr/num_legitimates;%false reject rate
        FAR(i)=100*fa/num_attackers;%false accept rate
    end 
end

%-------------- calculation of EER value
[~,tmp]=min(abs(FRR-FAR));
EER=FAR(tmp);tmpEER=tmp;


%-------------- calculation of the OP value
tmp2=find(OPvalue-FRR<=0);
tmpOP=min(length(tmp2)+1,length(FRR));
OP=FRR(tmpOP);far=FAR(tmpOP);frr=FRR(tmpOP);thre=x(tmpOP);
adr = 100 - far; % attack detection rate  = 100 - false accept rate

end

