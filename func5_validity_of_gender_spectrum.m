clear;clc;clf;close;

load('data\AAL2_FC.mat','net_HCP')
load('data\gender_continuum.mat','gc')  % gender continuum was obtained from func4_CV_and_testing_of_SVM_model function.

%% categorize all subjects into 3 groups according to their brain gender
for a = 1:11
    for b = 1:11
        n(a,b,1) = mean(net_HCP(a,b,gs<0.35));
        n(a,b,3) = mean(net_HCP(a,b,gs>0.65));
        n(a,b,2) = mean(net_HCP(a,b,gs>=0.35&gs<=0.65));
    end
end

%%
i_con = 0;              % edge counting index used to calculate color value
n_mono = 0;n_turn = 0;  % counting the number of edges that change monotonously
for a = 1:11
    for b = a:11
        hold on
        plot([1:3],squeeze(n(a,b,:)),'Marker','o','MarkerSize',4,'MarkerFaceColor',[(0+3.8*i_con)/255,(128+((-1)^i_con)*1.9*i_con)/255,1-(0+3.8*i_con)/255])
        i_con = i_con+1;
        if sign(n(a,b,2)-n(a,b,1))*sign(n(a,b,3)-n(a,b,2))>0
            n_mono = n_mono+1;
        else 
            n_turn = n_turn+1;
            fprintf('%.0f,%.0f\n',a,b) % print the edge that does not change monotonously
        end
    end
end
% ploting results
set(gca,'xtick',[1:3],'xticklabel',{});
xlim([0.8,3.2])

% network statistical analysis
i1 = find(gs<0.35);i2 = find(gs>=0.35 & gs<=0.65);i3 = find(gs>0.65);
group(i1) = 1;group(i2) = 2;group(i3) = 3;  % assigning group label in order to conduct ANOVA1 test

% edge wise one-way ANOVA test
for a = 1:11
    for b = a:11
        [p(a,b),tbl,~] = anova1(net_HCP(a,b,:),group,'off');
        F(a,b) = tbl{2,5};
    end
end

% multiple comparison adjustment
i = find(p~=0);p_vec = p(i);
sum(mafdr(p_vec,'BHFDR',1)<0.05);