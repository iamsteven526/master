
maxyy = 0;
for i = 1:10000000
    a1 = rand;
    b1 = rand;
    a2 = sqrt(1-a1*a1);
    b2 = sqrt(1-b1*b1);
    now = min([2-2*a1*b1-2*a2*b2,2-2*b1,3+2*a1-2*b1-2*a1*b1-2*a2*b2,3-2*a1-2*b1+2*a1*b1+2*a2*b2,2-2*a1]);
    if now > maxyy
        maxyy = now;
        nowa1 = a1;
        nowa2 = a2;
        nowb1 = b1;
        nowb2 = b2;
    end
end