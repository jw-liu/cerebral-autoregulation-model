function J_n = compute_outer_convergence(Q_n, Q_prev, A_n, A_prev)

N = size(Q_n, 1); J_sum = 0;
for w = 1:N
    dQ = max(abs(Q_n(w,:) - Q_prev(w,:))); dA = max(abs(A_n(w,:) - A_prev(w,:)));
    Q_ref = max(abs(Q_prev(w,:)));  A_ref = max(abs(A_prev(w,:)));
    if Q_ref > 0, J_sum = J_sum + dQ / Q_ref; end
    if A_ref > 0, J_sum = J_sum + dA / A_ref; end
end
J_n = J_sum / N;

end
