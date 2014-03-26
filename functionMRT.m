function wMRT = functionMRT(H,D)
%Calculates the maximum ratio transmission (MRT) beamforming vectors for a
%scenario where all or a subset of antennas transmit to each user.
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
%H  = Kr x Kt*Nt matrix with row index for users and column index
%     transmit antennas
%D  = Kt*Nt x Kt*Nt x Kr diagonal matrix. Element (j,j,k) is one if j:th
%     transmit antenna can transmit to user k and zero otherwise
%
%OUTPUT:
%wMRT = Kt*Nt x Kr matrix with normalized MRT beamforming



%Number of users
Kr = size(H,1);

%Total number of antennas
N = size(H,2);

%If D matrix is not provided, all antennas can transmit to everyone
if nargin<2
    D = repmat( eye(N), [1 1 Kr]);
end

%Pre-allocation of MRT beamforming
wMRT = zeros(size(H'));

%Computation of MRT, based on Definition 3.2
for k = 1:Kr
    channelvector = (H(k,:)*D(:,:,k))'; %Useful channel
    wMRT(:,k) = channelvector/norm(channelvector); %Normalization of useful channel
end
