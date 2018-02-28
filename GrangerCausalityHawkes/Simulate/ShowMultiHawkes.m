function ShowMultiHawkes(Seq, para, pattern)
% SHOWHAWKESM Plot intensities and event times of a Hawkes process
%             (Constant-unconditional-intensity first-order exponential-
%             decay M-variate Hawkes process; selected component series)
% Inputs      m   - indexes of selected component series, k*1 vector
%             t   - end time (Zero starting time assumed)
%             H   - process history, M*1 cell array; vector H{n} stores 
%                   event-occurrence times for component series n
%             par - process-parameters structure containing fields
%                   'mu'    - M*1 vector (unconditional intensities)
%                   'alpha' - M*M matrix (multiply the exponent terms)
%                   'beta'  - M*M matrix (degrees of the exponents)
% Example     See HAWKESDEMO
% Author      Dimitri Shvorob, dimitri.shvorob@vanderbilt.edu, 12/12/07


for n = 1:para.N
    seq = Seq{n};
    if ~isempty(seq)
        %x = linspace(0,para.T,100);


        set(gcf,'Color','w');
        lambda = zeros(para.U, size(seq,2));
        lambda(:,1)=para.mu;
        for i = 2:size(seq,2)
            lambda(:,i) = IntensityMultiHawkes(seq(1,i), seq, para, pattern);
        end

        figure
        for u = 1:para.U
            subplot(1,para.U,u)

            y = 0.05*max(lambda(u,:));
            e = find(seq(2,:)==u);
            if ~isempty(e)
                for j = e
                    line([seq(1,j) seq(1,j)],[0 y],'Color','k','LineWidth',2), hold on
                end    
                plot(seq(1,:),lambda(u,:),'Color',[0 .6 .6],'LineWidth',2)
                %axis([0 para.T 0 1.1*max(lambda(i,e))])
                ylabel('Intensity, \lambda(t)')
                xlabel(['Event-occurrence time (' num2str(length(e)) ' events total)'])
            end
        end
    end
end