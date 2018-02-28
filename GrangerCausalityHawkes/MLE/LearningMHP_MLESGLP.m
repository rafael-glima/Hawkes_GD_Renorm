function [Aest, muest, Landmark] = LearningMHP_MLESGLP( Seq, para, ...
                                                            alg, configure )

%
% Learning Nonparametric Hawkes Process Model with Sparse Group Lasso
%
% configure = [whether sparse, whether group-lasso, whether pairwise sim]
%

% initial 
Aest = 0.05*rand(alg.M, para.U, para.U);
muest = rand(para.U, 1)./para.U;


Landmark = (alg.T/alg.M)*(0:alg.M-1);
sigma = alg.sigma*ones(size(Landmark));




tic;
for o = 1:alg.outer
    
    
    for n = 1:alg.inner
        % initial objective function
        LLL = 0;
        
        Bmu = zeros(para.U,1);
        %BA = zeros(size(Aest));
        BmatA = alg.alphaS*ones(size(Aest));
        CmatA = zeros(size(Aest));
        AmatA = CmatA;
        Subgradient = CmatA;
        
        
        for u1 = 1:para.U
            for u2 = 1:para.U 
                tmp = Aest(:,u2,u1);
                if norm(tmp)==0
                    AmatA(:,u2,u1)=0;
                else
                    AmatA(:,u2,u1) = alg.alphaG/norm(tmp);
                end
            end
        end
        
        % E-step: evaluate the responsibility using the current parameters
        if configure(1)==1
            LLL = LLL + alg.alphaS*sum(abs(Aest(:)));
        end
        
        if configure(2)==1
            tmp = sqrt( sum( Aest.^2, 1 ) );
            LLL = LLL + alg.alphaG*sum(tmp(:));
        end
        
        if configure(3)==1
            valueP = 0;
            for np = 1:size(para.pair,1)
                ui = para.pair(np,1);
                uj = para.pair(np,2);
                vrow = Aest(:,:,ui) - Aest(:,:,uj);
                vcol = Aest(:,ui,:) - Aest(:,uj,:);
                
                Subgradient(:,:,ui) = Subgradient(:,:,ui) + 4*alg.alphaP*vrow;
                Subgradient(:,ui,:) = Subgradient(:,ui,:) + 4*alg.alphaP*vcol;
                
                AmatA(:,:,ui) = AmatA(:,:,ui) + 2*alg.alphaP;
                BmatA(:,:,ui) = BmatA(:,:,ui) - 2*alg.alphaP*Aest(:,:,uj);
                
                AmatA(:,:,uj) = AmatA(:,:,uj) + 2*alg.alphaP;
                BmatA(:,:,uj) = BmatA(:,:,uj) - 2*alg.alphaP*Aest(:,:,ui);
                
                AmatA(:,ui,:) = AmatA(:,ui,:) + 2*alg.alphaP;
                BmatA(:,ui,:) = BmatA(:,ui,:) - 2*alg.alphaP*Aest(:,uj,:);
                
                AmatA(:,uj,:) = AmatA(:,uj,:) + 2*alg.alphaP;
                BmatA(:,uj,:) = BmatA(:,uj,:) - 2*alg.alphaP*Aest(:,ui,:);
                
                valueP = valueP + 2*alg.alphaP*(sum(vrow(:).^2)+sum(vcol(:).^2));
            end
            LLL = LLL + valueP;    
        end
            
        for c = 1:length(Seq)
            seq = Seq{c};
            indt = seq(1,:)<para.T;
            seq = seq(:,indt);
            Time = seq(1,:);
            Event = seq(2,:);

            GK = Gkernel( para.T, Time, Landmark, sigma );
            
            
            Nc = length(Time);
            
            
            
            for i = 1:Nc
                
                

                ui = Event(i);
                
                BmatA(:,ui,:) = BmatA(:,ui,:)+...
                    double(Aest(:,ui,:)>0).*repmat( GK(:,i), [1,1,para.U] );
                Subgradient(:,ui,:) = Subgradient(:,ui,:)+...
                    double(Aest(:,ui,:)>0).*repmat( GK(:,i), [1,1,para.U] );
                ti = Time(i);             
                    
                lambdai = muest(ui);
                pii = muest(ui);
                pij = [];
                          
                SUM = pii;
                    
                if i>1
                    
                    tj = Time(1:i-1);
                    uj = Event(1:i-1);

                    DT = ti - tj;
                    index = find(DT<alg.T);
                    
                    if ~isempty(index)
                        tj = tj(index);
                        uj = uj(index);

                        gij = gkernel(ti, tj, Landmark, sigma);
                        auiuj = Aest(:, uj, ui);
                        pij = auiuj .* gij;

                        SUM = SUM+sum(pij(:));


                        lambdai = lambdai + sum(pij(:));
                    end
                 
                end

                LLL = LLL - log(lambdai);
                

                pii = pii./SUM;
                if i>1
                    pij = pij./SUM;
                    if ~isempty(pij) && sum(pij(:))>0
                        for j = 1:length(uj)
                            uuj = uj(j);
                            CmatA(:,uuj,ui) = CmatA(:,uuj,ui)-pij(:,j);                            
                        end
                    end
                end
                
                
                Bmu(ui) = Bmu(ui) + pii;
                
                
            end
            
            LLL = LLL + para.T.*sum(muest);
            LLL = LLL + sum( sum( GK.*sum(Aest(:,Event,:),3) ) );

            
            
        end
        
        Subgradient = Subgradient + (Aest>0).*CmatA./(Aest+eps);
        Subgradient(isnan(Subgradient)) = 0;
        Subgradient(isinf(Subgradient)) = 0;
        
        % M-step: update parameters
        mu = Bmu./(length(Seq)*para.T);
        if (alg.alphaG == 0 && alg.alphaP==0) || ...
                (configure(2)==0 && configure(3)==0)
            A = -CmatA./BmatA;
        else
            
            A = ( -BmatA + sqrt(BmatA.^2 - 4*AmatA.*CmatA) )./(2*AmatA);
            A(isnan(A))=0;
            A(isinf(A))=0;
        end
        
        
        
        
        % check convergence
        Err=sum(sum(sum(abs(A-Aest))))/(para.U*para.U*alg.M);
        fprintf('Outer=%d, Inner=%d, Obj=%f, RelErr=%f, Time=%0.2fsec\n',...
                o, n, LLL, Err, toc);
            
        if Err<alg.thres || (o==alg.outer && n==alg.inner)
            Aest = A;
            muest = mu;

            break;
        else
            Aest = A;
            muest = mu;

        end
        
    
    end

    if o<=alg.outer
        if configure(2)==1
            for u=1:para.U
                for v=1:para.U
                    if sum(abs(Aest(:,u,v)))>0
                        tmp = Aest(:,u,v) - alg.hardthres*Subgradient(:,u,v);

                        x = abs(tmp) - alg.hardthres*alg.alphaS;
                        x = sign(tmp).*double(x>=0).*x;

                        if norm(x) <= alg.hardthres*alg.alphaG
                            Aest(:,u,v)=0;
                        else
                            if configure(1)==1
                                res = 1-alg.hardthres*alg.alphaG/norm(x);
                                Aest(:,u,v) = double(res>=0).*abs(x);
                            end
                        end
                    end
                end
            end


        else
            for u=1:para.U
                for v=1:para.U
                    if sum(abs(Aest(:,u,v)))>0
                        if configure(1)==1
                            tmp = Aest(:,u,v) -...
                                alg.hardthres*Subgradient(:,u,v);
                            x = abs(tmp) - alg.hardthres*alg.alphaS;
                            x = sign(tmp).*double(x>=0).*x;
                            res = 1-alg.hardthres*alg.alphaG/norm(x);
                            Aest(:,u,v) = double(res>=0).*abs(x);
                        end
                    end
                end
            end
        end
    end
    Aest(Aest<0)=0;   
end
        

                
end


function v = gkernel( t, time, landmark, sigma )
    dt = t-time;
    dis = repmat(dt, [length(landmark), 1]) - repmat(landmark(:), [1,length(time)]);
    Sigma = repmat(sigma(:), [1,length(time)]);
    v = 1./(sqrt(2*pi).*Sigma) .* exp(-(dis.^2)./(2*Sigma.^2));
end

function v = Gkernel( T, time, landmark, sigma )
    dt = T-time;
    dis1 = repmat(dt, [length(landmark), 1]) - repmat(landmark(:), [1,length(time)]);
    dis2 = repmat(landmark(:), [1,length(time)]);
    Sigma = repmat(sigma(:), [1,length(time)]);
    v = 0.5*erf( dis1./(sqrt(2).*Sigma) )+0.5*erf( dis2./(sqrt(2).*Sigma) );
end
