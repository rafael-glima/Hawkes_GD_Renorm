%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate a set of event sequences based on multi-dimensional
% Hawkes processes
%
% Please cite our paper if you use our code
%
% Hongteng Xu, Mehrdad Farajtabar, and Hongyuan Zha. 
% "Learning granger causality for hawkes processes".
% International Conference on Machine Learning (ICML), 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% REMEMBER TO ADJUST MU AND N !!!!!!!!!!!!!!!!!

addpath('./Simulate/');
% time interval of sequence
para.N = 20;%5000;
para.T = 5000;
para.U = 1;
% intrinsic intensity matrix
%para.mu = ones(para.U,1); 
para.mu = .5 %rand(para.U, 1)./para.U; 
% distribution of components
para.pair = [1,2; 2,3; 4,5];
%para.freq = 0.1+blkdiag(0.2*ones(3), 0.1*ones(2));
para.freq = 0.1+blkdiag(0.2*ones(3), 0.2*ones(2));
%para.freq = zeros(5);
para.shift = double(para.freq<0.3);
%para.weight = blkdiag(1.0*ones(3), 1.0*ones(2)); 
para.weight = blkdiag(0.06*ones(3), 0.06*ones(2));
para.weight(1:3,4) = 0.06; %0.02;
para.weight(4,1:3) = 0.06; %0.02;
para.decayr = 0.2;
para.p = 11.;
decayr = para.decayr;

% figure
% hold on
% for i=1:para.U
%     for j=1:para.U
%         if para.shift(i,j)==0
%             dt=1/para.freq(i,j);
%             t=0:0.01:dt;        
%             plot(t, para.weight(i,j)*2*round(0.5*(1-cos(2*pi*para.freq(i,j)*t))));
%         else
%             dt=0.5/para.freq(i,j);
%             t=0:0.01:dt;        
%             plot(t, para.weight(i,j)*2*round(0.5*(1-cos(2*pi*para.freq(i,j)*t+pi*para.shift(i,j)))));
%         end
%     end
% end
% hold off

tic

Seq1 = SimMultiHawkes( para, 'exponential' );
Seq2 = SimMultiHawkes( para, 'powerlaw' );
Seq3 = SimMultiHawkes( para, 'q-exponential' );
Seq4 = SimMultiHawkes( para, 'rayleigh' );

 
save('4Kern_Renorm_20seq_T5000.mat','Seq1','Seq2','Seq3','Seq4','para'); 
%save('Sine_10seqT100000_oldfreq0.3.mat','Seq4','para');

%load Sine_10seqT100000_oldfreq0.3.mat
load 4Kern_Renorm_20seq_T5000.mat

figure
% subplot(2,2,1);
ShowKernel(para, 'exponential');
ShowMultiHawkes(Seq1, para, 'exponential');
%subplot(2,2,2);
ShowKernel(para, 'powerlaw');
ShowMultiHawkes(Seq2, para, 'powerlaw');
%subplot(2,2,3);
ShowKernel(para, 'q-exponential');
ShowMultiHawkes(Seq3, para, 'q-exponential');
% subplot(2,2,4);
ShowKernel(para, 'rayleigh');
ShowMultiHawkes(Seq4, para, 'rayleigh');

time = toc
