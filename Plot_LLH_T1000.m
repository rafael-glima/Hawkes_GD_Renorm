c = categorical({'EXP','PWL','QEXP','RAY'})

ax1 = subplot(1,4,1)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [-933.35 -915.78 -917.01 -917.337;-874.01 -869.30 -873.28 -874.01; -607.02 -576.74 -583.42 -584.68;-905.70 -905.70 -905.70 -905.70]
b = bar(c,score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b(1).FaceColor = [0 .5 .5]
b(2).FaceColor = [.5 0 .5]
% b(3).FaceColor = [0 .5 .5]
% b(4).FaceColor = [.5 0 .5]
% b(5).FaceColor = [0 .5 .5]
% b(6).FaceColor = [.5 0 .5]
% b(7).FaceColor = [0 .5 .5]
% b(8).FaceColor = [.5 0 .5]
b(1).EdgeColor = [0 .6 .6]
b(2).EdgeColor = [.6 0 .6]
% b(3).EdgeColor = [0 .6 .6]
% b(4).EdgeColor = [.6 0 .6]
% b(5).EdgeColor = [0 .6 .6]
% b(6).EdgeColor = [.6 0 .6]
% b(7).EdgeColor = [0 .6 .6]
% b(8).EdgeColor = [.6 0 .6]
%b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['EXP  (', char(949), ' = 0.1)'])
l = cell(1,4);
l{1}='MLE'; l{2}=['RF-MLE (', char(949), ' = 0.1)']; l{3}=['RF-MLE (', char(949), ' = 0.01)']; l{4}=['RF-MLE (', char(949), ' = 0.001)'];
legend(b,l)

ax2 = subplot(1,4,2)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [-1222.95 -1140.89 -1203.14 -1207.36; -1662.69 -1550.44 -1600.95 -1604.40;-979.52 -701.06 -740.48 -745.00;-1338.68 -1221.07 -1222.72 -1226.69]
b = bar(c,score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b(1).FaceColor = [0 .5 .5]
b(2).FaceColor = [.5 0 .5]
% b(3).FaceColor = [0 .5 .5]
% b(4).FaceColor = [.5 0 .5]
% b(5).FaceColor = [0 .5 .5]
% b(6).FaceColor = [.5 0 .5]
% b(7).FaceColor = [0 .5 .5]
% b(8).FaceColor = [.5 0 .5]
b(1).EdgeColor = [0 .6 .6]
b(2).EdgeColor = [.6 0 .6]
% b(3).EdgeColor = [0 .6 .6]
% b(4).EdgeColor = [.6 0 .6]
% b(5).EdgeColor = [0 .6 .6]
% b(6).EdgeColor = [.6 0 .6]
% b(7).EdgeColor = [0 .6 .6]
% b(8).EdgeColor = [.6 0 .6]
%b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['PWL  (', char(949), ' = 0.1)'])

ax3 = subplot(1,4,3)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [-781.88 -781.88 -781.88 -823.90;-848.25 -848.25 -848.25 -559.70; -779.39 -760.67 -760.67 -779.39;-850.53 -823.90 -781.88 -821.09]
b = bar(c,score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b(1).FaceColor = [0 .5 .5]
b(2).FaceColor = [.5 0 .5]
% b(3).FaceColor = [0 .5 .5]
% b(4).FaceColor = [.5 0 .5]
% b(5).FaceColor = [0 .5 .5]
% b(6).FaceColor = [.5 0 .5]
% b(7).FaceColor = [0 .5 .5]
% b(8).FaceColor = [.5 0 .5]
b(1).EdgeColor = [0 .6 .6]
b(2).EdgeColor = [.6 0 .6]
% b(3).EdgeColor = [0 .6 .6]
% b(4).EdgeColor = [.6 0 .6]
% b(5).EdgeColor = [0 .6 .6]
% b(6).EdgeColor = [.6 0 .6]
% b(7).EdgeColor = [0 .6 .6]
% b(8).EdgeColor = [.6 0 .6]
%b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['QEXP  (', char(949), ' = 0.1)'])

ax4 = subplot(1,4,4)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [-943.67 -934.11 -943.67 -943.67; -877.46 -877.46 -877.46 -877.46; -583.69 -566.26 -581.10 -583.69; -913.76 -913.76 -913.76 -913.76]
b = bar(c,score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b(1).FaceColor = [0 .5 .5]
b(2).FaceColor = [.5 0 .5]
% b(3).FaceColor = [0 .5 .5]
% b(4).FaceColor = [.5 0 .5]
% b(5).FaceColor = [0 .5 .5]
% b(6).FaceColor = [.5 0 .5]
% b(7).FaceColor = [0 .5 .5]
% b(8).FaceColor = [.5 0 .5]
b(1).EdgeColor = [0 .6 .6]
b(2).EdgeColor = [.6 0 .6]
%b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['RAY  (', char(949), ' = 0.1)'])

