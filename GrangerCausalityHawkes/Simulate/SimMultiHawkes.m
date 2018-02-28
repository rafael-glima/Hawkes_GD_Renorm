function Seq = SimMultiHawkes( para, pattern )

% Simulate event sequences of a mixture of multi-dimensional Hawkes process

Seq = cell(para.N, 1);

%tic
parfor n = 1:para.N


    n
    t=0;
    History = [];

    mt = SupIntensity(t, History, para, pattern);

    while t<para.T

        s = random('exp', 1/mt);
        U = rand;
        
        lambda_ts = IntensityMultiHawkes(t+s, History, para, pattern);
        mts = sum(lambda_ts);

        %fprintf('s=%f, v=%f\n', s, mts/mt);
        
        if t+s>para.T || U>mts/mt

            t = t+s;

        else

            u = rand*mts;
            sumIs = 0;
            for d=1:length(lambda_ts)
                sumIs = sumIs + lambda_ts(d);
                if sumIs >= u
                    break;
                end
            end
            index = d;


            t = t+s;
            History=[History,[t;index(1)]];

            

        end
        
        mt = SupIntensity(t, History, para, pattern);
        
    end
    Seq{n} = History;

    if mod(n, 10)==0 || n==para.N
%         fprintf('#seq=%d/%d, #event=%d, time=%.2fsec\n', ...
%             n, para.N, size(History,2));%, toc);
    end
    
        
end
    
    
