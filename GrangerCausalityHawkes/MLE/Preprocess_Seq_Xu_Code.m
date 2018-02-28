clear

N = 5000

addpath('./MLE');
addpath('./Simulate');

load 4Kern_1seqT100000_trunc.mat

Seq_3 = Seq3{1,1}

%Seq_3 = Seq_3(1:floor(Seq_3/N)*N)

Seq_3 = reshape(Seq_3,[N,length(Seq_3)/N])


