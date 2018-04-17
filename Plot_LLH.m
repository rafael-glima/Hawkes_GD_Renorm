c = categorical({'EXP','PWL','QEXP','RAY'})

ax1 = subplot(1,4,1)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [-933.35 -917.01 ;-874.01 -873.28; -607.02 -583.42 ;-905.70 -905.70]
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
l = cell(1,2);
l{1}='MLE'; l{2}='RF-MLE';    
legend(b,l)

ax2 = subplot(1,4,2)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [-1222.95 -1140.89; -1662.69 -1550.44 ;-979.52 -701.06 ;-1338.68 -1221.07]
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
score1 = [-781.88 -781.88 ;-848.25 -848.25; -779.39 -760.67 ;-850.53 -823.90]
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
score1 = [-943.67 -934.11; -877.46 -877.46; -583.69 -566.26; -913.76 -913.76]
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

