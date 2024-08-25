% ----------------------------------------------------------------------------
% This script solves a convex optimization problem to estimate the relative
% positions of mobile nodes and anchor nodes.
% We have m anchor nodes and n mobile nodes.
% Only partial distance information is available between the anchors, 
% as well as between the anchors and the mobile nodes. No distance 
% information is provided between the mobile nodes.

% Input:
% -------------------------------------------------------------------------
% Weight_E  = A matrix of zeros and ones with dimensions m*m. If the
% (i,j)-th entry of Weight_E is 1, it indicates that the distance 
% between anchor i and anchor j is observed.
% Weight_F = A matrix of zeros and ones with dimensions m*n. If the
% (i,j)-th entry of Weight_F is 1, it indicates that the distance 
% between anchor i and mobile node j is observed.
% -------------------------------------------------------------------------
% Output:
% -------------------------------------------------------------------------
% G_estimated = The estimated Gram matrix. One can then employ r-truncated
% eigendecomposition on this to obtain estimate of the configuration 
% of the points

% For a detailed explanation of the algorithm, refer to our paper below:
% Citation:
% Samuel Lichtenberg and Abiy Tasissa, “Localization from Structured Distance Matrices 
% via Low-Rank Matrix Recovery”, IEEE Transactions on Information Theory, 2024.
% arXiv link: https://arxiv.org/pdf/2311.18076
% ----------------------------------------------------------------------------
function G_estimated = AB_completion(D,Weight_E, Weight_F,m,n)
    % p
    p = m+n;
    % Construct the different blocks of the squared distance matrix
    E = D(1:m,1:m);
    F =  D(1:m,m+1:m+n);
    D_EF = [E F];
    % Construct the indices relevant to B and F
    [I,J] = find(Weight_F==1); 
    sz = [m, n]; 
    M = m*ones(length(I),1);
    ind_IJ = sub2ind(sz,I,J)+m*m;
    ind_mJ = sub2ind(sz,M,J)+m*m;
    % (I,I), (I,J) and (J,I) indices corresponding to observed entries
    % in the A-block
    szA = [m,m];
    ind_II = sub2ind(szA,I,I) ;
    [IA,JA] = find(Weight_E==1); 
    ind_IIA = sub2ind(szA,IA,IA);
    ind_JJA = sub2ind(szA,JA,JA);
    ind_IJA = sub2ind(szA,IA,JA);
    ind_JIA = sub2ind(szA,JA,IA);
    % CVX code
    cvx_begin 
    cvx_solver mosek
    cvx_precision high
    variable H(m,p)
    row_sum_H_A = sum(H(1:m,1:m),2); 
    tmp1 = m*H(ind_II)+trace(H(1:m,1:m))*ones(length(I),1)-2*row_sum_H_A(I);
    tmp2 = (0.5/m)*tmp1;
    tmp3 = -(0.5/m)*sum(E(m,:))*ones(length(I),1);
    minimize norm_nuc(H)
    H(1:m,1:m)==semidefinite(m,m);
    H(ind_IJ)-H(ind_mJ)-tmp2 == -0.5*(D_EF(ind_IJ)-D_EF(ind_mJ))+tmp3;
    H(ind_IIA)+H(ind_JJA)-2*H(ind_IJA) == E(ind_IJA);
    H(ind_IIA)+H(ind_JJA)-2*H(ind_JIA) == E(ind_JIA);
    sum(H)==zeros(1,p);
    cvx_end
    % Construct the estimated Gram matrix
    A_num = full(H(1:m,1:m));
    B_num = full(H(1:m,m+1:p));
    G_estimated = [A_num B_num; B_num' B_num'*(pinv(A_num,0.1))*B_num];
end