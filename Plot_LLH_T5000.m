c = categorical({'EXP','PWL','QEXP','RAY'})

ax1 = subplot(1,4,1)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [ -4676.61 -4599.52 ;-4409.35 -4386.74; -3133.68 -2998.176;-4479.72 -4479.72 ]
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
score1 = [-6427.87 -6049.23 ; -8732.92 -8324.04; -5162.15 -3550.52 ; -6198.44 -6031.68]
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
score1 = [-4263.04 -4764.31 ;-2246.73 -3028.74;-4834.63 ;-1898.07 ]
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
score1 = [;;;]
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

