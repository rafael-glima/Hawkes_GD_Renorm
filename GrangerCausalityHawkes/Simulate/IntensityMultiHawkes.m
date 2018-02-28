function lambda = IntensityMultiHawkes(t, History, para, pattern) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INTHAWKESM Intensity lambda(t) of a time-varying multi-dim Hawkes process
%             
% Inputs     t   - current interval 
%            History   - historical event sequences
%            T   - time interval
%            mu - U*1 vector (unconditional intensities)
%            w - parameters of decay function
%            Period, Shift - U*U parameters of infectivity matrix
% Outputs    lambda   - intensity 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(History)
    lambda = para.mu;
else
    Time = History(1, :);
    index = Time<=t;
    Time = Time(index);
    Event = History(2, index);

    lambda = para.mu;
    decayr = para.decayr;
    p = para.p;

    dt = t-Time;
    for u = 1:para.U
        for ui = 1:para.U
            if para.shift(u,ui)==0
                ind = dt<1/para.freq(u,ui) & Event == ui;
                lambda(u) = lambda(u) + sum( KernelFunc(dt(ind), para.weight(u, ui),...
                    para.freq(u, ui), para.shift(u, ui), pattern, decayr, p) );
            else
                ind = dt<0.5/para.freq(u,ui) & Event == ui;
                lambda(u) = lambda(u) + sum( KernelFunc(dt(ind), para.weight(u, ui),...
                    para.freq(u, ui), para.shift(u, ui), pattern, decayr, p) );
            end
        end
    

    end   
end

end

