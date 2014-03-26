%This Matlab script can be used to generate Figure 3 (a) and (b) in the
%lecture note:
%
%Emil Björnson, Mats Bengtsson, Björn Ottersten, “Optimal Multi-User
%Transmit Beamforming: Difficult Problem with a Simple Solution Structure,”
%IEEE Signal Processing Magazine, To appear.
%
%Download article: Address will be added here.
%
%The implementation of the heuristic beamforming schemes (MRT, ZFBF,
%transmit MMSE/regularized ZFBF/SLNR-MAX beamforming, and the corresponding
%power allocation is borrowed from the Matlab package of the following
%book:
%
%Emil Björnson, Eduard Jorswieck, “Optimal Resource Allocation in
%Coordinated Multi-Cell Systems,” Foundations and Trends in Communications
%and Information Theory, vol. 9, no. 2-3, pp. 113-381, 2013.
%
%The implementation of the branch-reduce-and-bound (BRB) algorithm, which
%computes the optimal beamforming, also originates from that book. The
%implemenation of this algorithm utilizes and requires CVX: http://cvxr.com/
%
%This is version 1.0 (Last edited: 2014-03-26)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the lecture note.


close all;
clear all;


%%Simulation parameters

rng('shuffle'); %Initiate the random number generators with a random seed
%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));

Nantennas = [4 12]; %Number of transmit antennas
K = 4; %Number of users

%Should the optimal beamforming be calculated? (true /false)
%(Note that it is the calculation of optimal beamforming is the main source
%of computational complexity. The running time with optimal beamforming
%takes hours or days, while it only takes a few minutes when only the
%heuristic beamforming schemes are computed. This shows clearly the need
%for simple heuristic beamforming!
computeOptimalBeamforming = true;

%Number of realizations in the Monte Carlo simulations
nbrOfMonteCarloRealizations = 2; %100;

%Combined channel matrix will be (K x K*N). This matrix gives the
%normalized variance of each channel element
channelVariances = [1 1 1 1];

%User weights for (unweighted) sum rate computation
weights = [1 1 1 1]'; ones(K,1);

%Range of SNR values
PdB = -10:1:30; %dB scale
P = 10.^(PdB/10); %Linear scale



%%System parameters - Special for optimal beamforming (using BRB algorithm)

%Optimal beamforming have its own (sparser) SNR range, since it is the
%Main source of computational complexity
PdB_BRB = -10:5:30; %SNR range (in dB) for optimal beamforming
P_BRB = 10.^(PdB_BRB/10); %SNR range for optimal beamforming

problemMode = 1; %Tells the BRB algorithm to maximize the (weighted) sum rate.

epsilonBRB = 0.01; %Accuracy of the optimal sum rate in the BRB algorithm
deltaBRB = 0.1; %Accuracy of the line searches in the BRB algorithm
maxIterations = 2000; %Maximal number of iterations of the algorithm
maxFuncEvaluations = 3000; %Maximal number of convex problems solved in the algorithm



%%Pre-allocation of matrices

%Matrices for saving sum rates with different beamforming stategies
sumRateMRT = zeros(length(P),nbrOfMonteCarloRealizations,length(Nantennas));
sumRateZFBF = zeros(length(P),nbrOfMonteCarloRealizations,length(Nantennas));
sumRateMMSE = zeros(length(P),nbrOfMonteCarloRealizations,length(Nantennas));

if computeOptimalBeamforming == true
    sumrateOPTIMAL = zeros(length(P_BRB),nbrOfMonteCarloRealizations,length(Nantennas));
end


%Go through the different number of transmit antennas
for n = 1:length(Nantennas)
    
    %Extract the current number of antennas
    N = Nantennas(n);
    
    %Pre-generation of Rayleigh fading channel realizations (unit variance)
    Hall = (randn(K,N,nbrOfMonteCarloRealizations)+1i*randn(K,N,nbrOfMonteCarloRealizations))/sqrt(2);
    
    %Go through all channel realizations
    for m = 1:nbrOfMonteCarloRealizations
        
    %Output the progress
    disp(['Progress: N = ' num2str(N) ', ' num2str(m) ' out of ' num2str(nbrOfMonteCarloRealizations) ' realizations.']);
        
        
        %Generate channel matrix for m:th realization
        H = repmat(sqrt(channelVariances)',[1 N]) .* Hall(:,:,m);
        
        
        %Compute normalized beamforming vectors for MRT
        wMRT = functionMRT(H);
        
        %Compute normalized beamforming vectors for ZFBF
        wZFBF = functionZFBF(H);
        
        
        %Pre-allocate matrix for saving good feasible starting points for the
        %BRB algorithm, based on heuristic beamforming
        if computeOptimalBeamforming == true
            bestfeasibleRates = zeros(K,length(P_BRB));
            pind_BRB = 1;
        end
        
        
        %Go through all transmit powers
        for pind = 1:length(P)
            
            %Compute normalized beamforming vectors for transmit MMSE
            %beamforming (which is the same as regularized ZFBF and SLNR-MAX
            %beamforming). Note that it varies with the transmit power.
            wSLNRMAX = functionSLNRMAX(H,P(pind)*ones(K,1));
            
            
            %Calculate power allocation with MRT (using Theorem 3.5 in [7])
            rhos = diag(abs(H*wMRT).^2)';
            powerAllocationMRT = functionHeuristicPowerAllocation(rhos,P(pind),weights);
            
            %Calculate sum rate with MRT
            W = kron(sqrt(powerAllocationMRT),ones(N,1)).*wMRT;
            channelGains = abs(H*W).^2;
            signalGains = diag(channelGains);
            interferenceGains = sum(channelGains,2)-signalGains;
            rates = log2(1+signalGains./(interferenceGains+1));
            sumRateMRT(pind,m,n) = weights'*rates;
            
            
            %Calculate power allocation with ZFBF (using Theorem 3.5 in [7])
            rhos = diag(abs(H*wZFBF).^2)';
            powerAllocationwZFBF = functionHeuristicPowerAllocation(rhos,P(pind),weights);
            
            %Calculate sum rate with ZFBF
            W = kron(sqrt(powerAllocationwZFBF),ones(N,1)).*wZFBF;
            channelGains = abs(H*W).^2;
            signalGains = diag(channelGains);
            interferenceGains = sum(channelGains,2)-signalGains;
            rates = log2(1+signalGains./(interferenceGains+1));
            sumRateZFBF(pind,m,n) = weights'*rates;
            
            
            %Calculate power allocation with transmit MMSE beamforming (using Theorem 3.5 in [7])
            rhos = diag(abs(H*wSLNRMAX).^2)';
            powerAllocationwSLNRMAX_sumrate = functionHeuristicPowerAllocation(rhos,P(pind),weights);
            
            %Calculate sum rate with transmit MMSE beamforming
            W = kron(sqrt(powerAllocationwSLNRMAX_sumrate),ones(N,1)).*wSLNRMAX;
            channelGains = abs(H*W).^2;
            signalGains = diag(channelGains);
            interferenceGains = sum(channelGains,2)-signalGains;
            rates = log2(1+signalGains./(interferenceGains+1));
            sumRateMMSE(pind,m,n) = weights'*rates;
            
            
            
            %Save the rates of transmit MMSE beamforming to use as starting
            %point when calculating Optimal beamforming
            if computeOptimalBeamforming == true && P(pind)==P_BRB(pind_BRB)
                bestfeasibleRates(:,pind_BRB) = rates;
                pind_BRB = pind_BRB+1;
            end
            
        end
        
        if computeOptimalBeamforming == true
            
            for pind = 1:length(P_BRB)
                %[pind length(P_BRB)]  %Simple output to keep track of progress
                
                %Definition of a total power constraint
                L = 1;
                
                %The matrix-square root of the weighting matrix
                Qsqrt = sqrt(1/P_BRB(pind))*eye(N);
                
                %Normalized limit of the total transmit power
                q = ones(L,1);
                
                
                %The BRB algorithm searches in a box with lower corner in the
                %origin. The upper corner is the utopia point, where each user
                %gets the rate that it would get if it was the only user in the
                %system. This point is computed below.
                origin = zeros(K,1);
                
                %Computation of the utopia point using MRT, which is optimal
                %when there is only one active user
                utopia = zeros(K,1);
                for k = 1:K
                    utopia(k) = log2(1+abs(H(k,:)*wMRT(:,k))^2*P_BRB(pind));
                end
                
                %This matrix is similar to diagonal matrices defined in the
                %section "Multiple Cooperating Base Stations" and says that all
                %antennas transmit to all users.
                D = repmat(eye(N),[1 1 K]);
                
                %Obtain the beamforming that maximizes the sum rate using the
                %BRB algorithm from [7].
                bestFeasibleBRB = functionBRBalgorithm_cvx(H,D,Qsqrt,q,origin,utopia,weights,deltaBRB,epsilonBRB,maxIterations,maxFuncEvaluations,bestfeasibleRates(:,pind),problemMode);
                
                %Save the performance of the optimal beamforming
                sumrateOPTIMAL(pind,m,n) = weights'*bestFeasibleBRB;
            end
            
        end
        
    end
    
end




%Plot simulation results
for n = 1:length(Nantennas)

figure; hold on; box on;

if computeOptimalBeamforming==true
    plot(PdB_BRB,mean(sumrateOPTIMAL(:,:,n),2),'k:','LineWidth',1);
end

plot(PdB,mean(sumRateMMSE(:,:,n),2),'r','LineWidth',1);
plot(PdB,mean(sumRateZFBF(:,:,n),2),'b--','LineWidth',1);
plot(PdB,mean(sumRateMRT(:,:,n),2),'k-.','LineWidth',1);

if computeOptimalBeamforming == true
    legend('Optimal Beamforming','Transmit MMSE/Regul. ZFBF','ZFBF','MRT','Location','NorthWest');
elseif computeOptimalBeamforming == false
    legend('Transmit MMSE/Regul. ZFBF','ZFBF','MRT','Location','NorthWest');
end

title(['N = ' num2str(Nantennas(n)) ' antennas']);

xlabel('Average SNR [dB]')
ylabel('Average Sum Rate [bit/channel use]');

end