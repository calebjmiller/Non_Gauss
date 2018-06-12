function cov_mat = real_get_cov_Gaussian_needlet(beta, b_mat,b_mat_deriv, Npix, A_1,A_2)

% number of basis functions
m = size(b_mat, 2)-1;
eta = [0 beta(1:m)];
std_vec = exp(b_mat*reshape(eta, m+1, 1)); 
std_vec_deriv = b_mat_deriv*reshape(eta, m+1, 1).*std_vec;
sigma_sq = beta(m+1:end);

[N, M] = size(A_1);

DA = zeros(N, M);
for i = 1:N
    DA(i, :) = std_vec(i)*A_1(i, :) - std_vec_deriv(i)*A_2(i,:);
end

len_j = length(Npix);
st = zeros(len_j, 1);
en = zeros(len_j, 1);
for j = 1:len_j
    st(j) = sum(Npix(1:j))-Npix(j)+1;
    en(j) = sum(Npix(1:j));
end

cov_mat = zeros(N);
for j = 1:len_j
    range = st(j):en(j);
    DA_sub = DA(:, range);
    cov_mat = cov_mat+sigma_sq(j)*(DA_sub*DA_sub');
end

end