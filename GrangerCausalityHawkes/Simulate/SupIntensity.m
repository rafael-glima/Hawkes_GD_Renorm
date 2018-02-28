function mt = SupIntensity(t, History, para, pattern) 

%
% Compute the super bound of intensity function 
%



if isempty(History)
    mt = sum(para.mu);
else
    Time = History(1, :);
    index = Time<=t;
    Time = Time(index);
    Event = History(2, index);

    tmp = para.freq;
%     if strcmp(pattern,'exponential') || strcmp(pattern,'powerlaw')
%         dtmax = Time;
%     else
%         dtmax = round(1/min(tmp(:)));
%     end
    dtmax = round(1/min(tmp(:)));
    M=50;

    
    MT = para.mu*ones(1, M);
    for m=1:M
        dt = t+(m-1)*dtmax/M - Time;
        
        for u = 1:para.U
            for ui = 1:para.U
                if para.shift(u,ui)==0
                    ind =  dt<1/para.freq(u,ui) & Event == ui ;
                    %if ~isempty(ind)
                    MT(u,m) = MT(u,m) + sum( KernelFunc(dt(ind), ...
                        para.weight(u,ui), para.freq(u,ui), para.shift(u,ui), pattern, para.decayr, para.p) );
                    %end
                else
                    ind =  dt<0.5/para.freq(u,ui) & Event == ui ;
                    %if ~isempty(ind)
                    MT(u,m) = MT(u,m) + sum( KernelFunc(dt(ind), ...
                        para.weight(u,ui), para.freq(u,ui), para.shift(u,ui), pattern, para.decayr, para.p) );
                    %end
                end
            end
        end
        
    end  
    mt = max(sum(MT));
end



end




