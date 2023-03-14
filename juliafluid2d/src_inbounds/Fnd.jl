module Fnd
export MP5

@inline function minmod(x,y)
    sgn = sign(x)
    return  sgn*max(min(abs(x),sgn*y),0.0)
end

@inline function MUSCL3(qm,q,
    qp,q2p)
    k = 1/3
    b = (3-k)/(1-k)
    dp  = qp - q
    dm  = q - qm
    ddp = minmod(dp, b * dm)
    ddm = minmod(dm, b * dp)
    qL  = q + 0.25 * ((1.0 - k) * ddm + (1.0 + k) * ddp)
    dp  = q2p - qp
    dm  = qp - q
    ddp = minmod(dp, b * dm)
    ddm = minmod(dm, b * dp)
    qR = qp - 0.25 * ((1.0 - k) * ddp + (1.0 + k) * ddm)
    return qL,qR
end

@inline function minmod4(
    x1,x2,
    x3,x4
    )
    return 0.5 * (sign(x1) + sign(x2)) * abs(
        0.5 * (sign(x1) + sign(x3)) *
        0.5 * (sign(x1) + sign(x4))) *
        min(abs(x1),abs(x2),abs(x3),abs(x4) );
end

@inline function median(x,y,z)
    return x + minmod(y - x, z - x)
end

@inline function MP5_sub(q2m,qm,
        q,qp,
        q2p
        )
    alp = 2.0
    qL = (2.0 * q2m - 13.0 * qm + 47.0 * q + 27.0 * qp - 3.0 * q2p) / 60.0
    qMP = q + minmod(qp - q, alp * (q - qm))

    if (qL-q)*(qL-qMP)>1.0e-10
        qLL = qL
        dm = q2m + q - 2.0 * qm
        d = qm + qp - 2.0 * q;
        dp = q + q2p - 2.0 * qp;
        dMm = minmod4(4.0 * dm - d, 4.0 * d - dm, dm, d);
        dMp = minmod4(4.0 * d - dp, 4.0 * dp - d, d, dp);
        qUL = q + alp * (q - qm);
        qMD = 0.5 * (q + qp) - 0.5 * dMp;
        qLC = q + 0.5 * (q - qm) + 4.0 / 3.0 * dMm;
        qmin = max(min(q,qp,qMD ), min(q,qUL,qLC ));
        qmax = min(max( q,qp,qMD ), max(q,qUL,qLC ));
        qL = median(qLL, qmin, qmax)

    end

    return qL
end


function MP5(q2m,qm,
    q,qp,
    q2p,q3p,
    )

    qL = MP5_sub(q2m, qm, q, qp, q2p)
    qR = MP5_sub(q3p, q2p, qp, q, qm)

    return qL,qR

end

end