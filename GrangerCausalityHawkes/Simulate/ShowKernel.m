function ShowKernel(para, pattern)

decayr = para.decayr;
p = para.p;
figure



for i=1:para.U
    for j=1:para.U
        subplot(para.U, para.U, para.U*(i-1)+j)
        if para.shift(i,j)==0
            dt=1/para.freq(i,j);
            dt
            %axis([-1.0 3.5 0.0 0.2])
            t=0:0.01:dt;             
            plot(t, KernelFunc( t, para.weight(i,j), para.freq(i,j), para.shift(i,j), pattern, decayr, p),'-k','LineWidth',5);
            %axis off
        else
            dt=0.5/para.freq(i,j);
            dt
            %axis([-1.0 3.5 0.0 0.2])
            t=0:0.01:dt;        
            plot(t, KernelFunc( t, para.weight(i,j), para.freq(i,j), para.shift(i,j), pattern, decayr, p),'-k','LineWidth',5);
            %axis off
        end
    end
end
