function wSLNRMAX = functionSLNRMAX(H,eta,D)
%Calculates the Signal-to-leakage-and-noise ratio maximizing (SLNR-MAX)
%beamforming for a scenario where all or a subset of antennas transmit 
%to each user. Note that SLNR-MAX is also known as regularized zero-forcing
%beamforming and transmit MMSE beamforming
%
%The references to definitions and equations refer to the following book:
%
%Emil Björnson, Eduard Jorswieck, “Optimal Resource Allocation in
%Coordinated Multi-Cell Systems,” Foundations and Trends in Communications
%and Information Theory, vol. 9, no. 2-3, pp. 113-381, 2013.
%
%This is version 1.1. (Last edited: 2014-03-26)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%INPUT:
%H   = Kr x Kt*Nt matrix with row index for users and column index
%      transmit antennas
%eta = Kr x 1 vector with SNR^(-1) like parameter of this user 
%D   = Kt*Nt x Kt*Nt x Kr diagonal matrix. Element (j,j,k) is one if j:th
%      transmit antenna can transmit to user k and zero otherwise
%
%OUTPUT:
%wSLNRMAX = Kt*Nt x Kr matrix with normalized SLNR-MAX beamforming



%Number of users
Kr = size(H,1);

%Total number of antennas
N = size(H,2);

%If eta vector is not provided, all values are set to unity
if nargin<2
    eta = ones(Kr,1);
end

%If D matrix is not provided, all antennas can transmit to everyone
if nargin<3
    D = repmat( eye(N), [1 1 Kr]);
end

%Pre-allocation of MRT beamforming
wSLNRMAX = zeros(size(H'));

%Computation of SLNR-MAX, based on Definition 3.5
for k = 1:Kr
    effectivechannel = (H*D(:,:,k))'; %Effective channels
    projectedchannel = (eye(N)/eta(k)+effectivechannel*effectivechannel')\effectivechannel(:,k); %Compute zero-forcing based on channel inversion
    wSLNRMAX(:,k) = projectedchannel/norm(projectedchannel);  %Normalization of zero-forcing direction
end
