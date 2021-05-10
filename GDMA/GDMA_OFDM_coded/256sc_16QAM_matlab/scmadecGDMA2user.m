function LLR = scmadecGDMA2user(y, h, N0)

M = 2;
V = 2;

% Factor graph calculation

N   = size(y, 2);
LLR = zeros(log2(M)*V, N);

Noise = 1/N0;

for jj = 1:N

    % Step 1: Initial calculations
    f = zeros(M, M);

     % non-zero elements, paths
    for m1 = 1:M
        for m2 = 1:M
            d = y(jj)-(2*(1.5-m1)*h(1,jj)+2*(1.5-m2)*h(2,jj));
            f(m1,m2) = exp(-Noise*sum(real(d)^2+imag(d)^2));
        end
    end

% Step 3: LLR calculation
    Q = zeros(M, V);
    Q(1,1) = sum(sum(f(1,:)));
    Q(2,1) = sum(sum(f(2,:)));
    Q(1,2) = sum(sum(f(:,1)));
    Q(2,2) = sum(sum(f(:,2)));    

    LLR_tmp = zeros(log2(M)*V, 1);
    for k = 1:V
        LLR_tmp(k) = log(Q(1,k)/Q(2,k));
    end
    LLR(:,jj) = LLR_tmp;
end
