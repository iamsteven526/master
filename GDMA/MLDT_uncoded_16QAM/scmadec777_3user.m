function LLR = scmadec777_3user(y, h, N0)

M = 16;
V = 3;

% Factor graph calculation

N   = size(y, 2);
LLR = zeros(V, log2(M)*N);

Noise = 1/N0;

QAM_point = qammod([0:15],16,'UnitAveragePower',true,'InputType','integer')';

for jj = 1:N

    % Step 1: Initial calculations
    f = zeros(M, M, M);
    
     % non-zero elements, paths
    for m1 = 1:M
        for m2 = 1:M
            for m3 = 1:M
                d = y(jj)-(QAM_point(m1)*h(1,jj)+QAM_point(m2)*h(2,jj)+QAM_point(m3)*h(3,jj));
                f(m1,m2,m3) = exp(-Noise*sum(real(d)^2+imag(d)^2));
            end
        end
    end

% Step 3: LLR calculation
    Q = zeros(M, V);
    for mm = 1:M
        Q(mm,1) = sum(sum(f(mm,:,:)));
        Q(mm,2) = sum(sum(f(:,mm,:)));
        Q(mm,3) = sum(sum(f(:,:,mm)));
    end
    %Q(1,1) = sum(sum(f(1,:,:)));
    %Q(2,1) = sum(sum(f(2,:,:)));
    %Q(1,2) = sum(sum(f(:,1,:)));
    %Q(2,2) = sum(sum(f(:,2,:)));    
    %Q(1,3) = sum(sum(f(:,:,1)));
    %Q(2,3) = sum(sum(f(:,:,2)));

    LLR_tmp = zeros(V, log2(M)*1); % temp variable for parallel work

    for k = 1:V
        LLR_tmp(k,1) = log(sum(Q(1:8,k))/sum(Q(9:16,k)));
        LLR_tmp(k,2) = log((sum(Q(1:4,k))+sum(Q(9:12,k)))/(sum(Q(5:8,k))+sum(Q(13:16,k))));
        LLR_tmp(k,3) = log((sum(Q(1:4:13,k))+sum(Q(2:4:14,k)))/(sum(Q(3:4:15,k))+sum(Q(4:4:16,k))));
        LLR_tmp(k,4) = log(sum(Q(1:2:15,k))/sum(Q(2:2:16,k)));
    end
    LLR(:,4*(jj-1)+1:4*jj) = LLR_tmp;
end
