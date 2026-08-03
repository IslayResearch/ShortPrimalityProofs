/* short.gp -- compute a short ECPP certificate for a probable prime p >= 5.
 * Written by Fable 5.  A toy implementation, in the style of oneshot.gp; the goal is
 * clarity, not speed.
 *
 * Format: github.com/AndrewVSutherland/ShortPrimalityProofs.  The certificate is the flat
 * sequence  (p_0, A_0, x_0, o_0, A_1, x_1, o_1, ..., A_k, x_k, o_k)  with o_i = m_i*p_{i+1},
 * p_{k+1} = 1 and n = ceil(log_2 p_0) FIXED for the whole chain: on E_{A_i} : y^2 = x^3 +
 * A_i x^2 + x over F_{p_i} the point with x-coordinate x_i has order exactly o_i, where m_i
 * is n^2-smooth, p_{i+1} is a prime in (n^2, sqrt(p_i)) proven prime by the next level, and
 * L_i < o_i < r_i L_i with L_i = (p_i^{1/4}+1)^2 and r_i the least prime divisor of m_i.
 *
 * Method.  For each level we search random curves E_A/F_{p_i}, compute N = #E_A by SEA
 * (ellcard), and look for a divisor o | N inside the narrow window (L_i, r_i L_i) -- note
 * o is only about sqrt(p_i), so most of N is discarded.  A usable o must be of the shape
 * (n^2-smooth) * (one prime), so we trial-divide N up to n^2 to get its smooth part s,
 * and then look for the single large prime among the factors of the rough part N/s.  That
 * last step is the expensive one: the rough part is factored under an `alarm` time budget
 * (SC_tlim seconds) and the curve is abandoned if the budget runs out -- an early abort
 * that keeps the search on curves whose order factors easily.  A level with a fully
 * n^2-smooth o (p_{i+1} = 1) ends the chain.
 *
 * Usage:
 *     echo 'printshort(nextprime(10^30))' | gp -q short.gp
 *     echo 'SC_tlim = 60; printshort(nextprime(10^100))' | gp -q short.gp
 */

default(parisizemax, 2^30);                              \\ allow the stack to grow for large SEA computations
default(factor_add_primes, 1);                           \\ retain factors found inside timed-out factorint calls

SC_curves = 0;                                           \\ curve orders tried (all levels)
SC_seacalls = 0;                                         \\ expensive ellcard calls (one call gives two twists)
SC_tlim = 20;                                            \\ seconds allowed per rough-part factorization
                                                         \\ (SC_tlim = 0 disables PARI factorint; saved
                                                         \\  prefactoring or external ECM may still run)
SC_factorflags = 0;                                      \\ factorint flags; 9 is a fast partial-factor mode
SC_factorrounds = 3;                                     \\ continue only after a factor pass shrinks its residual
SC_deepfactorbits = 0;                                   \\ use a larger budget after reducing R to this size
SC_deepfactortlim = 0;                                   \\ larger factor budget (0 disables adaptive factoring)
SC_prefactortlim = 0;                                    \\ optional saved partial-factor pass (0 disables)
SC_prefactorflags = 8;                                   \\ factorint flags for the saved partial pass
SC_pm1b1 = 0;                                            \\ optional saved GMP P-1 stage-1 bound
SC_pp1b1 = 0;                                            \\ optional saved GMP P+1 stage-1 bound
SC_ecmb1 = 0;                                            \\ optional GMP-ECM stage-1 bound (0 disables)
SC_ecmcurves = 0;                                        \\ GMP-ECM curves tried per rough cofactor
SC_ecmrounds = 1;                                        \\ saved ECM factors extracted per rough cofactor
SC_msievebits = 0;                                       \\ optional msieve residual threshold (0 disables)
SC_msieveseconds = 0;                                    \\ wall-time limit for each msieve call
SC_workerkey = 0;                                        \\ unique runner key for temporary external files
SC_curvefamily = 0;                                      \\ 0 random; 1/2/3 known 4/8/16-torsion;
                                                         \\ 4 known 3-torsion; 5 X_1(27) with point 4
SC_rootsmoothbits = 0;                                   \\ root smooth-part threshold (0 disables)
SC_candidatebits = 0;                                    \\ save root orders with at least this many smooth bits
SC_candidatefile = "";                                   \\ append-only resumable root-order log
SC_levelfile = "";                                       \\ append-only proven root-level log
SC_rootp = 0;
SC_branchcurves = 64;                                    \\ backtrack after this many curves below a descent
                                                         \\ (0 retains the original unlimited commitment)
SC_maxcurves = 0;                                        \\ 0 = unlimited; else give up after this many curves
SC_progress = 0;                                         \\ emit cumulative counters after every tested order
SC_seators = 0;                                          \\ ellsea torsion filter (negative filters the twist);
                                                         \\ 0 uses unfiltered ellcard
SC_seaaborts = 0;                                        \\ SEA calls rejected early by SC_seators
SC_factorattempts = 0;                                   \\ rough-part factorization calls
SC_factoraborts = 0;                                     \\ factorization calls that hit an alarm/error
SC_descents = 0;                                         \\ nonterminal levels found
SC_backtracks = 0;                                       \\ child subtrees abandoned
SC_ecmattempts = 0;                                      \\ rough cofactors sent to GMP-ECM
SC_ecmfactors = 0;                                       \\ GMP-ECM calls that found a factor
SC_qcandidates = 0;                                      \\ prime factors in the allowed q range
SC_viableqcandidates = 0;                                \\ q candidates with s*q above the lower window
SC_windowcandidates = 0;                                 \\ (m,q) pairs in the certificate window
SC_residualbits = 0;                                     \\ bit length entering the full factor pass
SC_maxsmoothbits = 0;                                    \\ largest smooth-part bit length observed
SC_smoothrejects = 0;                                    \\ root orders below SC_rootsmoothbits
SC_deepfactorattempts = 0;                               \\ factor calls receiving the adaptive larger budget
SC_pm1attempts = 0; SC_pm1factors = 0;                   \\ saved GMP P-1 calls/factors
SC_pp1attempts = 0; SC_pp1factors = 0;                   \\ saved GMP P+1 calls/factors
SC_msieveattempts = 0; SC_msievefactors = 0;             \\ bounded msieve calls that ran/found factors
SC_exhaustedorders = 0;                                  \\ fully split orders with no admissible factor
SC_factorrecoveries = 0;                                 \\ timed-out factorizations that still saved progress
SC_cmtests = 0; SC_cmlastD = 0;                          \\ discriminants tested/current |D|
SC_resumeattempts = 0;                                   \\ queued root checkpoints actually attempted
SC_smoothbound = 0;                                      \\ bound represented by cached prime product
SC_smoothprimeproduct = 1;                               \\ product of every prime through SC_smoothbound

scprogress() = {
  if(SC_progress,
    warning("progress curves=", SC_curves, " sea_calls=", SC_seacalls,
            " sea_aborts=", SC_seaaborts, " factor_attempts=", SC_factorattempts,
            " factor_aborts=", SC_factoraborts, " descents=", SC_descents,
            " backtracks=", SC_backtracks, " ecm_attempts=", SC_ecmattempts,
            " ecm_factors=", SC_ecmfactors, " q_candidates=", SC_qcandidates,
            " viable_q=", SC_viableqcandidates, " window_candidates=", SC_windowcandidates,
            " residual_bits=", SC_residualbits, " max_smooth_bits=", SC_maxsmoothbits,
            " smooth_rejects=", SC_smoothrejects,
            " deep_factor_attempts=", SC_deepfactorattempts,
            " pm1_attempts=", SC_pm1attempts, " pm1_factors=", SC_pm1factors,
            " pp1_attempts=", SC_pp1attempts, " pp1_factors=", SC_pp1factors,
            " msieve_attempts=", SC_msieveattempts,
            " msieve_factors=", SC_msievefactors,
            " exhausted_orders=", SC_exhaustedorders,
            " factor_recoveries=", SC_factorrecoveries,
            " cm_tests=", SC_cmtests, " cm_last_D=", SC_cmlastD,
            " resume_attempts=", SC_resumeattempts));
};

/* Prepare the square-free product of every prime through B.  Repeated gcds with
 * this product extract every prime power from an order exactly, while replacing
 * tens of thousands of single-prime divisions with a few quasi-linear big-int
 * operations.  The product is under 200 KiB even for the 300-digit target. */
scpreparesmooth(B) = {
  if(B == SC_smoothbound, return());
  SC_smoothprimeproduct = 1;
  forprime(q = 2, B, SC_smoothprimeproduct *= q);
  SC_smoothbound = B;
};

/* n^2-smooth part s of N and the rough cofactor r = N/s. */
smoothpart(N, B) = {
  my(s = 1, r = N, g);
  if(B == SC_smoothbound && SC_smoothprimeproduct > 1,
    while(1,
      g = gcd(r, SC_smoothprimeproduct);
      if(g == 1, break);
      s *= g; r /= g
    );
    return([s, r])
  );
  forprime(q = 2, B, while(r % q == 0, r /= q; s *= q));
  [s, r];
};

/* PARI can retain prime factors discovered before alarm() expires.  Re-running a
 * bounded trial factorization recovers those private-table entries in milliseconds,
 * so a timeout no longer discards successful ECM/Rho work.  Return [matrix,timed_out]. */
sctimedfactor(R, seconds, flags) = {
  my(F = alarm(seconds, factorint(R, flags)), timedout = 0);
  if(type(F) == "t_ERROR",
    timedout = 1;
    \\ factor_add_primes already placed factors found before the alarm in PARI's
    \\ add-prime table.  A limit of 2 tests only that table (R is already n^2-rough),
    \\ avoiding an accidental second trial division through 2^24.
    F = factor(R, 2)
  );
  [F, timedout];
};

scbound(p) = sqrtint(p) + 1 + sqrtint(4 * sqrtint(p));   \\ integer form of L = (p^{1/4}+1)^2

/* A small quadratic nonresidue modulo p. */
scnonsquare(p) = {
  my(d);
  if(p % 4 == 3, return(p-1));
  d = 2; while(kronecker(d, p) != -1, d++);
  d;
};

/* Construct a Montgomery coefficient from Sutherland's optimized X_1(27)
 * model.  The finite-field model/map implementation is adapted from
 * IslayResearch/OneShotSEA (MIT; see THIRD_PARTY_NOTICES.md), which pins the
 * published formula at https://math.mit.edu/~drew/X1/X1opt27new.txt.
 * We retain only models with full rational 2-torsion and a rational
 * point of order 4.  The distinguished order-27 point then makes the selected
 * curve order divisible by 216; the paired quadratic twist is still tested for
 * the price of the same SEA call. */
scx127A(p) = {
  my(u, v, um1, u2pu1, u5, u6, c0, c1, c2, c3, c4, c5, c6, P, vs,
     g, x, x2, x3, g2, g3, g4, yn, yd, y, xy, x2y, rden, sden, r, s,
     rs, r2s, a1, a2, a3, b2, b4, b6, aa, bb, roots, ei, ej, ek,
     point4, deriv, d, A, X = 'X);
  while(1,
    u = Mod(random(p), p); if(!u, next);
    um1 = u - 1; u2pu1 = u^2 + u + 1; u5 = u^5; u6 = u^6;
    c6 = um1^2;
    c5 = um1^2 * (u^3 + 2);
    c4 = -um1^2 * (u^5 + 2*u^4 - 2*u^3 - u^2 - 2*u - 1);
    c3 = u*um1 * (u^6 - 3*u^5 - 4*u^4 + u^3 + u^2 + 3*u - 2);
    c2 = u*um1*u2pu1 * (3*u^4 - 4*u^3 - 2*u^2 + u - 1);
    c1 = 3*u5*um1*u2pu1;
    c0 = u6*u2pu1;
    P = lift(c0 + c1*X + c2*X^2 + c3*X^3 + c4*X^4 + c5*X^5 + c6*X^6);
    vs = iferr(polrootsmod(P, p), e, []);
    for(vi = 1, #vs,
      v = Mod(vs[vi], p);
      if(!u*(v+1), next);
      g = -1/u; x = v/(u*(v+1));
      x2 = x^2; x3 = x^3; g2 = g^2; g3 = g^3; g4 = g^4;
      yn = g4*x + g4 + g3*x - 2*g2*x2 - g*x3 + g + x;
      yd = g4 + g3*x - g2*x2 - g*x2 + g + x;
      if(!yd, next);
      y = yn/yd; xy = x*y; x2y = x2*y;
      rden = x2y-x; sden = xy;
      if(!rden || !sden, next);
      r = (x2y-xy+y-1)/rden;
      s = (xy-y+1)/sden;
      rs = r*s; r2s = r*rs;
      a1 = s-rs+1; a2 = rs-r2s; a3 = a2;
      b2 = a1^2 + 4*a2; b4 = a1*a3; b6 = a3^2;
      c4 = b2^2 - 24*b4;
      c6 = -b2^3 + 36*b2*b4 - 216*b6;
      aa = -27*c4; bb = -54*c6;
      if(4*aa^3 + 27*bb^2 == 0, next);
      roots = iferr(polrootsmod(lift(X^3 + aa*X + bb), p), e, []);
      if(#roots != 3, next);
      point4 = 0;
      for(i = 1, 3,
        ei = Mod(roots[i], p);
        ej = Mod(roots[if(i == 1, 2, 1)], p);
        ek = Mod(roots[if(i == 3, 2, 3)], p);
        if(kronecker(lift(ei-ej), p) == 1 &&
           kronecker(lift(ei-ek), p) == 1,
          point4 = 1; break)
      );
      if(!point4, next);
      for(i = 1, 3,
        ei = Mod(roots[i], p); deriv = 3*ei^2 + aa;
        if(kronecker(lift(deriv), p) != 1, next);
        d = sqrt(deriv); A = lift(3*ei/d);
        if((A^2-4) % p, return(A))
      )
    )
  )
};

/* Choose A for a random Montgomery curve, optionally with known torsion. */
sccurveA(p) = {
  my(u, x, t, A, rhs, P, roots, X = 'X);
  if(SC_curvefamily == 1,
    u = random(p);
    return((u^2 - 2) % p)                                \\ (1,u) has order 4
  );
  if(SC_curvefamily == 2,
    while(1,
      x = random(p); if(!x, next);
      A = lift(Mod((x^2 - 1)^2 - 4*x*(x^2 + 1), p) / Mod(4*x^2, p));
      rhs = (x * (x^2 + A*x + 1)) % p;
      if(kronecker(rhs, p) == 1, return(A))               \\ doubling x maps x to 1, then to 0
    )
  );
  if(SC_curvefamily == 3,
    while(1,
      t = random(p); if(!t, next);
      A = lift(Mod((t^2 - 1)^2 - 4*t*(t^2 + 1), p) / Mod(4*t^2, p));
      rhs = (t * (t^2 + A*t + 1)) % p;
      if(kronecker(rhs, p) != 1, next);                   \\ t is the x-coordinate of 8-torsion
      P = (X^2 - 1)^2 - 4*t*X*(X^2 + A*X + 1);
      roots = polrootsmod(P, p);
      for(i = 1, #roots,
        u = lift(roots[i]); if(!u, next);
        rhs = (u * (u^2 + A*u + 1)) % p;
        if(kronecker(rhs, p) == 1, return(A))             \\ u halves t, giving 16-torsion
      )
    )
  );
  if(SC_curvefamily == 4,
    while(1,
      t = random(p); if(!t, next);
      A = lift(Mod((t^2 - 1)^2 - 4*t^2*(t^2 + 1), p) / Mod(4*t^3, p));
      rhs = (t * (t^2 + A*t + 1)) % p;
      if(kronecker(rhs, p) == 1, return(A))               \\ doubling fixes x=t, giving 3-torsion
    )
  );
  if(SC_curvefamily == 5, return(scx127A(p)));
  random(p);
};

/* A point of order exactly o on E (N = #E, o | N, fo = the primes of o), or 0 if none found. */
scpoint(E, N, o, fo) = {
  my(Q, ok);
  for(t = 1, 64,
    Q = ellmul(E, random(E), N/o);                       \\ order(Q) divides o
    if(#Q == 1, next);                                   \\ Q = O, resample
    if(#ellmul(E, Q, o) != 1, next);                     \\ also safe for a proposed CM order
    ok = 1;
    for(i = 1, #fo, if(#ellmul(E, Q, o/fo[i]) == 1, ok = 0; break));
    if(ok, return(Q))
  );
  0;
};

/* Find an admissible (o,q) pair in an already factored rough part, without needing
 * a curve model.  This is used to screen cheap CM orders before constructing their
 * (potentially much more expensive) Hilbert class polynomial. */
scpickorders(n2, L, rt, s, dv, F) = {
  my(q, m, o, C = List());
  if(type(F) != "t_COL", return(0));
  for(j = 1, #F, q = F[j];
    if(q > n2 && q < rt && ispseudoprime(q),
      SC_qcandidates++;
      if(s*q > L, SC_viableqcandidates++);
      for(i = 1, #dv, m = dv[i]; o = m * q;
        if(m > 1 && o > L && o < factor(m)[1,1] * L,
          SC_windowcandidates++;
          listput(C, [o, q, concat(factor(m)[,1], [q]~)])))));
  if(#C, Vec(C), 0);
};

/* Test the pseudoprime entries in a factor list as the large prime q in a descent. */
sctryfactors(p, n2, A, xden, E, N, L, rt, s, dv, F) = {
  my(C = scpickorders(n2, L, rt, s, dv, F), c, Q);
  if(type(C) == "t_VEC",
    for(k = 1, #C, c = C[k];
    Q = scpoint(E, N, c[1], c[3]);
      if(Q != 0, return([A, lift(Q[1]/Mod(xden, p)), c[1], c[2]])))
  );
  0;
};

/* Append an exact resumable root order.  Repeated records for the same order checkpoint
 * progressively smaller composite residuals after saved factoring stages. */
sclogcandidate(p, A, xden, N, s, R) = {
  if(p == SC_rootp && SC_candidatebits && SC_candidatefile != "" &&
     #binary(s) >= SC_candidatebits,
    write(SC_candidatefile, Str(p, " ", A, " ", xden, " ", N, " ", s, " ", R))
  );
};

/* Save a CM order using xden=0 as an unambiguous record tag and |D| in the A field.
 * The Montgomery model is deliberately deferred until factoring exposes a usable q. */
sclogcmcandidate(p, D, N, s, R) = {
  if(p == SC_rootp && SC_candidatebits && SC_candidatefile != "" &&
     #binary(s) >= SC_candidatebits,
    write(SC_candidatefile, Str(p, " ", -D, " 0 ", N, " ", s, " ", R))
  );
};

/* Once the exact-order point test has succeeded, preserve the complete root level.
 * Five-field records are distinguished from the six-field order checkpoints above. */
scloglevel(p, lev) = {
  if(p == SC_rootp && SC_levelfile != "",
    write(SC_levelfile, Str(p, " ", lev[1], " ", lev[2], " ", lev[3], " ", lev[4]))
  );
};

/* Preserve CM orders that have passed factoring and the certificate-window screen before
 * beginning the potentially long class-polynomial reconstruction.  Seven fields distinguish
 * these records from six-field rough orders and five-field complete root levels. */
sclogcmscreen(p, D, N, C) = {
  if(p == SC_rootp && SC_levelfile != "" && type(C) == "t_VEC",
    for(i = 1, #C,
      write(SC_levelfile, Str(p, " ", -D, " ", N, " ", C[i][1], " ", C[i][2], " 0 0")))
  );
};

/* Mark each attempted CM screen as deterministically exhausted after the class polynomial
 * was constructed successfully but supplied no matching Montgomery level. */
sclogcmtombstone(p, D, N, C) = {
  if(p == SC_rootp && SC_levelfile != "" && type(C) == "t_VEC",
    for(i = 1, #C,
      write(SC_levelfile, Str(p, " ", -D, " ", N, " ", C[i][1], " ",
                              C[i][2], " 1 0")))
  );
};

/* Run optional external GMP-ECM curves and return its factor list. */
scecmfactors(R) = {
  my(out, F);
  if(!SC_ecmb1 || !SC_ecmcurves, return(0));
  out = iferr(externstr(Str("printf '%s\\n' ", R, " | ecm -q -c ", SC_ecmcurves,
                            " -one ", SC_ecmb1)), e, []);
  if(type(out) != "t_VEC" || #out < 1, return(0));
  F = apply(x -> eval(x), strsplit(out[1], " "));
  if(type(F) != "t_VEC" || #F < 1, return(0));
  F~;
};

/* Run one optional external GMP P-1 or P+1 pass and return its factor list. */
scspecialfactors(R, method, B1) = {
  my(out, F);
  if(!B1, return(0));
  out = iferr(externstr(Str("printf '%s\\n' ", R, " | ecm -q -", method,
                            " -one ", B1)), e, []);
  if(type(out) != "t_VEC" || #out < 1, return(0));
  F = apply(x -> eval(x), strsplit(out[1], " "));
  if(type(F) != "t_VEC" || #F < 1, return(0));
  F~;
};

/* Run a wall-time-bounded one-thread msieve pass and retain every emitted factor. */
scmsievefactors(R) = {
  my(out, fields, q, U = R, F = List(), base);
  if(!SC_msievebits || !SC_msieveseconds || #binary(R) > SC_msievebits, return(0));
  base = Str("/tmp/short-msieve-", SC_workerkey, "-", R % 1000000007);
  out = iferr(externstr(Str("rm -f ", base, ".dat ", base, ".log; (timeout ",
                            SC_msieveseconds, " msieve -q -t 1 -s ", base,
                            ".dat -l ", base, ".log ", R, " || true); rm -f ",
                            base, ".dat ", base, ".log")), e, []);
  if(type(out) != "t_VEC", return(0));
  for(i = 1, #out,
    fields = strsplit(out[i], ": ");
    if(#fields == 2,
      q = iferr(eval(fields[2]), e, 0);
      if(type(q) == "t_INT" && q > 1 && U % q == 0,
        listput(F, q);
        while(U % q == 0, U /= q)
      )
    )
  );
  if(U > 1, listput(F, U));
  Vec(F)~;
};

/* Try one curve order that has already been computed.  A is the certificate parameter and
 * xden maps the model's x-coordinate to the certificate coordinate x/xden. */
sctryorder(p, n2, A, xden, E, N, L, rt, {R0 = 0}) = {
  my(sr, s, R, oldR, dv, F, EF, PF, TF, U, lev, m, o, q, Q, factortlim);
    sr = smoothpart(N, n2); s = sr[1]; R = sr[2];
    if(R0,
      if(R % R0, error("short: saved residual does not divide the rough part"));
      R = R0
    );
    SC_maxsmoothbits = max(SC_maxsmoothbits, #binary(s));
    sclogcandidate(p, A, xden, N, s, R);
    if(p == SC_rootp && SC_rootsmoothbits && #binary(s) < SC_rootsmoothbits,
      SC_smoothrejects++;
      return(0)
    );
    dv = divisors(s);

    \\ (a) terminal level: o = m is n^2-smooth all by itself, and lands in the window
    for(i = 1, #dv, m = dv[i];
      if(m > L && m < factor(m)[1,1] * L,
        Q = scpoint(E, N, m, factor(m)[,1]);
        if(Q != 0, return([A, lift(Q[1]/Mod(xden, p)), m, 1]))));

    \\ (b) cheap descent: if the saved rough residual is itself prime, try it directly.
    \\     This also handles reduced resume checkpoints, where prime factors already tested and
    \\     discarded mean that s*R need not equal N.  A rejected prime residual exhausts the order.
    if(R > 1 && ispseudoprime(R),
      lev = sctryfactors(p, n2, A, xden, E, N, L, rt, s, dv, [R]~);
      if(lev != 0, return(lev));
      sclogcandidate(p, A, xden, N, s, 1);
      SC_exhaustedorders++;
      return(0)
    );

    \\ (c) general descent: o = m*q with q any prime factor of the rough part, m | s, m > 1.
    \\     Here the discarded cofactor N/o is unconstrained, which is far likelier -- but q must
    \\     be dug out of R, so we factor R under a time budget and abandon slow curves.
    if((SC_tlim > 0 || SC_deepfactortlim > 0 || SC_prefactortlim > 0 ||
        (SC_msievebits > 0 && SC_msieveseconds > 0) ||
        SC_pm1b1 > 0 || SC_pp1b1 > 0 ||
        (SC_ecmb1 > 0 && SC_ecmcurves > 0)) && R > 1 &&
        s * min(R \ (n2+1), rt) > L,
      if(SC_pm1b1 > 0,
        SC_pm1attempts++;
        EF = scspecialfactors(R, "pm1", SC_pm1b1);
        if(type(EF) == "t_COL",
          lev = sctryfactors(p, n2, A, xden, E, N, L, rt, s, dv, EF);
          if(lev != 0, return(lev));
          U = 1;
          for(j = 1, #EF, if(!ispseudoprime(EF[j]), U *= EF[j]));
          if(#EF > 1, SC_pm1factors++);
          if(U == 1,
            sclogcandidate(p, A, xden, N, s, 1);
            SC_exhaustedorders++; return(0));
          if(#EF > 1, R = U; sclogcandidate(p, A, xden, N, s, R))
        )
      );
      if(SC_pp1b1 > 0 && s * min(R \ (n2+1), rt) > L,
        SC_pp1attempts++;
        EF = scspecialfactors(R, "pp1", SC_pp1b1);
        if(type(EF) == "t_COL",
          lev = sctryfactors(p, n2, A, xden, E, N, L, rt, s, dv, EF);
          if(lev != 0, return(lev));
          U = 1;
          for(j = 1, #EF, if(!ispseudoprime(EF[j]), U *= EF[j]));
          if(#EF > 1, SC_pp1factors++);
          if(U == 1,
            sclogcandidate(p, A, xden, N, s, 1);
            SC_exhaustedorders++; return(0));
          if(#EF > 1, R = U; sclogcandidate(p, A, xden, N, s, R))
        )
      );
      if(SC_ecmb1 > 0 && SC_ecmcurves > 0 && s * min(R \ (n2+1), rt) > L,
        for(round = 1, SC_ecmrounds,
          if(s * min(R \ (n2+1), rt) <= L, break);
          SC_ecmattempts++;
          EF = scecmfactors(R);
          if(type(EF) != "t_COL", break);
          lev = sctryfactors(p, n2, A, xden, E, N, L, rt, s, dv, EF);
          if(lev != 0, return(lev));
          U = 1;
          for(j = 1, #EF, if(!ispseudoprime(EF[j]), U *= EF[j]));
          if(#EF > 1, SC_ecmfactors++);
          if(U == 1,
            sclogcandidate(p, A, xden, N, s, 1);
            SC_exhaustedorders++; return(0));
          \\ With `ecm -one`, a no-factor run prints R itself.  That parses as
          \\ a one-entry vector and must end the saved-factor loop; another
          \\ round is useful only after ECM actually shrank the residual.
          if(#EF <= 1, break);
          R = U;
          sclogcandidate(p, A, xden, N, s, R)
        )
      );
      if(SC_prefactortlim > 0 && s * min(R \ (n2+1), rt) > L,
        SC_factorattempts++;
        TF = sctimedfactor(R, SC_prefactortlim, SC_prefactorflags);
        PF = TF[1];
        if(TF[2], SC_factoraborts++);
        if(type(PF) == "t_MAT",
          lev = sctryfactors(p, n2, A, xden, E, N, L, rt, s, dv, PF[,1]);
          if(lev != 0, return(lev));
          U = 1;
          for(j = 1, matsize(PF)[1],
            if(!ispseudoprime(PF[j,1]), U *= PF[j,1]^PF[j,2]));
          if(TF[2] && U < R, SC_factorrecoveries++);
          if(U == 1,
            sclogcandidate(p, A, xden, N, s, 1);
            SC_exhaustedorders++; return(0));
          R = U;
          sclogcandidate(p, A, xden, N, s, R)
        )
      );
      SC_residualbits = #binary(R);
      if(SC_msievebits > 0 && SC_msieveseconds > 0 &&
         s * min(R \ (n2+1), rt) > L &&
         SC_residualbits <= SC_msievebits,
        SC_msieveattempts++;
        F = scmsievefactors(R);
        if(type(F) == "t_COL",
          lev = sctryfactors(p, n2, A, xden, E, N, L, rt, s, dv, F);
          if(lev != 0, return(lev));
          U = 1;
          for(j = 1, #F, if(!ispseudoprime(F[j]), U *= F[j]));
          if(#F > 1, SC_msievefactors++);
          if(U == 1,
            sclogcandidate(p, A, xden, N, s, 1);
            SC_exhaustedorders++; return(0));
          if(#F > 1, R = U; sclogcandidate(p, A, xden, N, s, R))
        );
        SC_residualbits = #binary(R)
      );
      factortlim = SC_tlim;
      if(SC_deepfactorbits > 0 && SC_deepfactortlim > factortlim &&
         SC_residualbits <= SC_deepfactorbits,
        factortlim = SC_deepfactortlim;
        SC_deepfactorattempts++
      );
      if(factortlim > 0 && s * min(R \ (n2+1), rt) > L,
        for(fr = 1, SC_factorrounds,
          if(s * min(R \ (n2+1), rt) <= L, break);
          oldR = R;
          SC_factorattempts++;
          TF = sctimedfactor(R, factortlim, SC_factorflags);
          PF = TF[1];
          if(TF[2], SC_factoraborts++);
          F = if(type(PF) == "t_MAT", PF[,1], 0);
                                                         \\ flag 9 may leave a composite cofactor, filtered below
          lev = sctryfactors(p, n2, A, xden, E, N, L, rt, s, dv, F);
          if(lev != 0, return(lev));
          if(type(F) != "t_COL", break);
          U = 1;
          for(j = 1, #F, if(!ispseudoprime(F[j]), U *= F[j]));
          if(TF[2] && U < oldR, SC_factorrecoveries++);
          if(U == 1,
            sclogcandidate(p, A, xden, N, s, 1);
            SC_exhaustedorders++; return(0));
          if(U >= oldR, break);
          R = U;
          SC_residualbits = #binary(R);
          sclogcandidate(p, A, xden, N, s, R)
        )));
    /* Every still-composite factor of R is larger than n2.  If even the
     * largest possible remaining prime, R/(n2+1), cannot clear the
     * certificate window with s, this order is permanently exhausted. */
    if(R == 1 || s * min(R \ (n2+1), rt) <= L,
      sclogcandidate(p, A, xden, N, s, 1);
      SC_exhaustedorders++
    );
  0;
};

/* One level: search random curves over F_p for [A, x0, o, p_next].  p_next = 1 ends the chain.
 * If d is a quadratic nonresidue, y^2 = x^3+d*A*x^2+d^2*x is the quadratic twist of E_A,
 * and its model coordinate X maps to certificate coordinate x=X/d on d*y^2=x^3+A*x^2+x.
 * The two orders sum to 2p+2, so test both for the cost of one SEA computation. */
sclevel(p, n2, {stopcurve = 0}) = {
  my(L = scbound(p), rt = sqrtint(p), d = scnonsquare(p), A, E, Et, N, lev);
  while(1,
    if(SC_maxcurves && SC_curves >= SC_maxcurves, return(0));
    if(stopcurve && SC_curves >= stopcurve, return(0));
    A = sccurveA(p); E = ellinit([0, A, 0, 1, 0], p);
    SC_curves++;
    if(#E == 0, next);                                   \\ singular (A = +-2 mod p)
    SC_seacalls++;
    N = if(SC_seators, ellsea(E, SC_seators), ellcard(E)); \\ SEA, optionally with early abort
    if(!N,
      SC_seaaborts++;
      scprogress();
      next
    );
    lev = sctryorder(p, n2, A, 1, E, N, L, rt);
    scprogress();
    if(lev != 0, return(lev));

    if(SC_maxcurves && SC_curves >= SC_maxcurves, return(0));
    if(stopcurve && SC_curves >= stopcurve, return(0));
    Et = ellinit([0, (d*A)%p, 0, (d^2)%p, 0], p);
    SC_curves++;
    lev = sctryorder(p, n2, A, d, Et, 2*p+2-N, L, rt);
    scprogress();
    if(lev != 0, return(lev))
  );
};

/* Screen an exact CM order before constructing its curve.  CM gives the two possible
 * orders p+1+-t essentially for free; only orders with a sufficiently large smooth part
 * pay the configured saved-factor portfolio.  Return every admissible
 * [o,q,prime divisors of o] candidate, or zero. */
scscreenorder(p, n2, N, L, rt, {D = 0}, {R0 = 0}) = {
  my(sr = smoothpart(N, n2), s = sr[1], R = sr[2], oldR, dv, m, F, EF, PF, TF, U, C,
     factortlim, terminal = List());
  /* Every certificate curve has the rational 2-torsion point (0,0), so its
   * full group order is even.  An odd CM order can never yield this model. */
  if(N % 2, return(0));
  if(R0,
    if(R % R0, error("short: saved CM residual does not divide the rough part"));
    R = R0
  );
  SC_maxsmoothbits = max(SC_maxsmoothbits, #binary(s));
  if(D, sclogcmcandidate(p, D, N, s, R));
  if(p == SC_rootp && SC_rootsmoothbits && #binary(s) < SC_rootsmoothbits,
    SC_smoothrejects++;
    return(0)
  );
  dv = divisors(s);
  for(i = 1, #dv, m = dv[i];
    if(m > L && m < factor(m)[1,1] * L,
      SC_windowcandidates++;
      listput(terminal, [m, 1, factor(m)[,1]])));
  if(#terminal, return(Vec(terminal)));
  if(R > 1 && ispseudoprime(R),
    C = scpickorders(n2, L, rt, s, dv, [R]~);
    if(type(C) == "t_VEC", return(C));
    if(D, sclogcmcandidate(p, D, N, s, 1));
    SC_exhaustedorders++;
    return(0)
  );
  if(SC_pm1b1 > 0 && R > 1 && s * min(R \ (n2+1), rt) > L,
    SC_pm1attempts++;
    EF = scspecialfactors(R, "pm1", SC_pm1b1);
    if(type(EF) == "t_COL",
      C = scpickorders(n2, L, rt, s, dv, EF);
      if(type(C) == "t_VEC", return(C));
      U = 1;
      for(j = 1, #EF, if(!ispseudoprime(EF[j]), U *= EF[j]));
      if(#EF > 1, SC_pm1factors++);
      if(U == 1,
        if(D, sclogcmcandidate(p, D, N, s, 1));
        SC_exhaustedorders++; return(0));
      if(#EF > 1,
        R = U;
        if(D, sclogcmcandidate(p, D, N, s, R)))
    )
  );
  if(SC_pp1b1 > 0 && R > 1 && s * min(R \ (n2+1), rt) > L,
    SC_pp1attempts++;
    EF = scspecialfactors(R, "pp1", SC_pp1b1);
    if(type(EF) == "t_COL",
      C = scpickorders(n2, L, rt, s, dv, EF);
      if(type(C) == "t_VEC", return(C));
      U = 1;
      for(j = 1, #EF, if(!ispseudoprime(EF[j]), U *= EF[j]));
      if(#EF > 1, SC_pp1factors++);
      if(U == 1,
        if(D, sclogcmcandidate(p, D, N, s, 1));
        SC_exhaustedorders++; return(0));
      if(#EF > 1,
        R = U;
        if(D, sclogcmcandidate(p, D, N, s, R)))
    )
  );
  if(SC_ecmb1 > 0 && SC_ecmcurves > 0 && R > 1 &&
     s * min(R \ (n2+1), rt) > L,
    for(round = 1, SC_ecmrounds,
      if(s * min(R \ (n2+1), rt) <= L, break);
      SC_ecmattempts++;
      EF = scecmfactors(R);
      if(type(EF) != "t_COL", break);
      C = scpickorders(n2, L, rt, s, dv, EF);
      if(type(C) == "t_VEC", return(C));
      U = 1;
      for(j = 1, #EF, if(!ispseudoprime(EF[j]), U *= EF[j]));
      if(#EF > 1, SC_ecmfactors++);
      if(U == 1,
        if(D, sclogcmcandidate(p, D, N, s, 1));
        SC_exhaustedorders++; return(0));
      if(#EF <= 1, break);
      R = U;
      if(D, sclogcmcandidate(p, D, N, s, R))
    )
  );
  if(SC_prefactortlim > 0 && R > 1 && s * min(R \ (n2+1), rt) > L,
    SC_factorattempts++;
    TF = sctimedfactor(R, SC_prefactortlim, SC_prefactorflags);
    PF = TF[1];
    if(TF[2], SC_factoraborts++);
    if(type(PF) == "t_MAT",
      C = scpickorders(n2, L, rt, s, dv, PF[,1]);
      if(type(C) == "t_VEC", return(C));
      U = 1;
      for(j = 1, matsize(PF)[1],
        if(!ispseudoprime(PF[j,1]), U *= PF[j,1]^PF[j,2]));
      if(TF[2] && U < R, SC_factorrecoveries++);
      if(U == 1,
        if(D, sclogcmcandidate(p, D, N, s, 1));
        SC_exhaustedorders++; return(0));
      R = U;
      if(D, sclogcmcandidate(p, D, N, s, R))
    )
  );
  SC_residualbits = #binary(R);
  if(SC_msievebits > 0 && SC_msieveseconds > 0 && R > 1 &&
     s * min(R \ (n2+1), rt) > L &&
     SC_residualbits <= SC_msievebits,
    SC_msieveattempts++;
    F = scmsievefactors(R);
    if(type(F) == "t_COL",
      C = scpickorders(n2, L, rt, s, dv, F);
      if(type(C) == "t_VEC", return(C));
      U = 1;
      for(j = 1, #F, if(!ispseudoprime(F[j]), U *= F[j]));
      if(#F > 1, SC_msievefactors++);
      if(U == 1,
        if(D, sclogcmcandidate(p, D, N, s, 1));
        SC_exhaustedorders++; return(0));
      if(#F > 1,
        R = U;
        if(D, sclogcmcandidate(p, D, N, s, R)))
    );
    SC_residualbits = #binary(R)
  );
  factortlim = SC_tlim;
  if(SC_deepfactorbits > 0 && SC_deepfactortlim > factortlim &&
     SC_residualbits <= SC_deepfactorbits,
    factortlim = SC_deepfactortlim;
    SC_deepfactorattempts++
  );
  if(factortlim > 0 && R > 1 && s * min(R \ (n2+1), rt) > L,
    for(fr = 1, SC_factorrounds,
      if(s * min(R \ (n2+1), rt) <= L, break);
      oldR = R;
      SC_factorattempts++;
      TF = sctimedfactor(R, factortlim, SC_factorflags);
      PF = TF[1];
      if(TF[2], SC_factoraborts++);
      F = if(type(PF) == "t_MAT", PF[,1], 0);
      C = scpickorders(n2, L, rt, s, dv, F);
      if(type(C) == "t_VEC", return(C));
      if(type(F) != "t_COL", break);
      U = 1;
      for(j = 1, #F, if(!ispseudoprime(F[j]), U *= F[j]));
      if(TF[2] && U < oldR, SC_factorrecoveries++);
      if(U == 1,
        if(D, sclogcmcandidate(p, D, N, s, 1));
        SC_exhaustedorders++; return(0));
      if(U >= oldR, break);
      R = U;
      SC_residualbits = #binary(R);
      if(D, sclogcmcandidate(p, D, N, s, R))
    )
  );
  /* Any unresolved factor is n2-rough.  Once R/(n2+1) is too small,
   * no refinement of this composite residual can expose a usable q. */
  if(R == 1 || s * min(R \ (n2+1), rt) <= L,
    if(D, sclogcmcandidate(p, D, N, s, 1));
    SC_exhaustedorders++
  );
  0;
};

/* Construct a Montgomery model with CM discriminant D and find the certified point.
 * For j(E_A)=256(A^2-3)^3/(A^2-4), solve first for A^2, then try both square roots.
 * The explicit [o]Q=O check in scpoint safely distinguishes the two proposed CM orders. */
scmontgomerylevel(p, D, N, C) = {
  my(H, ws, jinv, zs, z, a, A, E, Q, c,
     failed = 0, inv = if(D % 3, 5, 0), X = 'X);
  H = iferr(polclass(D, inv), e, 0);
  if(type(H) != "t_POL", return(-1));
  ws = iferr(polrootsmod(H, p), e, 0);
  if(type(ws) != "t_COL" && type(ws) != "t_VEC", return(-1));
  for(i = 1, #ws,
    jinv = if(inv == 5, ws[i]^3, ws[i]);
    zs = iferr(polrootsmod(256*(X-3)^3 - lift(jinv)*(X-4), p), e,
               failed = 1; []);
    if(type(zs) != "t_COL" && type(zs) != "t_VEC", next);
    for(j = 1, #zs,
      z = lift(zs[j]);
      if(kronecker(z, p) == -1, next);
      a = lift(sqrt(Mod(z, p)));
      forstep(sign = -1, 1, 2,
        A = lift(Mod(sign*a, p));
        E = ellinit([0, A, 0, 1, 0], p);
        if(#E == 0, next);
        for(k = 1, #C, c = C[k];
          Q = scpoint(E, N, c[1], c[3]);
          if(Q != 0, return([A, lift(Q[1]), c[1], c[2]])))
      )
    )
  );
  if(failed, -1, 0);
};

/* Recursively complete a chain, abandoning a child subtree at the absolute `stopcurve`.
 * A zero stop is unlimited (used only for the root), while every descent receives at most
 * the SC_branchcurves budget.  This avoids committing a worker forever to an unlucky child. */
scchain(p, n2, stopcurve) = {
  if(p == 1, return([]));
  my(lev, tail, childstop);
  while(1,
    if(SC_maxcurves && SC_curves >= SC_maxcurves, return(0));
    if(stopcurve && SC_curves >= stopcurve, return(0));
    lev = sclevel(p, n2, stopcurve);
    if(lev == 0, return(0));
    scloglevel(p, lev);
    if(lev[4] != 1, SC_descents++);
    childstop = if(SC_branchcurves, SC_curves + SC_branchcurves, 0);
    if(stopcurve && (!childstop || stopcurve < childstop), childstop = stopcurve);
    tail = scchain(lev[4], n2, childstop);
    if(type(tail) == "t_VEC", return(concat([lev[1], lev[2], lev[3]], tail)));
    SC_backtracks++
  );
};

/* Reset search counters before a fresh or resumed root search. */
screset(p) = {
  SC_rootp = p;
  scpreparesmooth((#binary(p))^2);
  SC_curves = 0; SC_seacalls = 0; SC_seaaborts = 0;
  SC_factorattempts = 0; SC_factoraborts = 0; SC_descents = 0; SC_backtracks = 0;
  SC_ecmattempts = 0; SC_ecmfactors = 0;
  SC_qcandidates = 0; SC_viableqcandidates = 0; SC_windowcandidates = 0;
  SC_residualbits = 0; SC_maxsmoothbits = 0;
  SC_smoothrejects = 0; SC_deepfactorattempts = 0;
  SC_pm1attempts = 0; SC_pm1factors = 0; SC_pp1attempts = 0; SC_pp1factors = 0;
  SC_msieveattempts = 0; SC_msievefactors = 0;
  SC_exhaustedorders = 0;
  SC_factorrecoveries = 0;
  SC_cmtests = 0; SC_cmlastD = 0;
  SC_resumeattempts = 0;
};

/* The full chain: returns the flat sequence (p, A_0, x_0, o_0, ..., A_k, x_k, o_k). */
shortcert(p, {reset = 1}) = {
  if(!ispseudoprime(p), error("short: p is composite"));
  if(p < 5, error("short: need p >= 5"));
  my(n = #binary(p), n2 = n^2, tail);
  if(reset, screset(p));
  tail = scchain(p, n2, 0);
  if(type(tail) != "t_VEC", return(0));
  concat([p], tail);
};

/* Try the two exact orders supplied by one CM trace, construct a Montgomery model only
 * after an order passes factoring, and recursively finish the certificate. */
sccmtrace(p, n2, L, rt, D, t) = {
  my(N, C, lev, tail, childstop);
  forstep(sign = -1, 1, 2,
    if(SC_maxcurves && SC_curves >= SC_maxcurves, return(0));
    N = p + 1 + sign*t;
    /* Montgomery certificate curves have rational 2-torsion.  Reject the
     * impossible odd order before doing any saved factoring work. */
    if(N % 2, next);
    SC_curves++;
    C = scscreenorder(p, n2, N, L, rt, D);
    scprogress();
    if(C == 0, next);
    sclogcmscreen(p, D, N, C);
    lev = scmontgomerylevel(p, D, N, C);
    if(type(lev) != "t_VEC",
      if(lev == 0, sclogcmtombstone(p, D, N, C));
      next);
    scloglevel(p, lev);
    if(lev[4] != 1, SC_descents++);
    childstop = if(SC_branchcurves, SC_curves + SC_branchcurves, 0);
    tail = scchain(lev[4], n2, childstop);
    if(type(tail) == "t_VEC",
      return(concat([p, lev[1], lev[2], lev[3]], tail)));
    SC_backtracks++
  );
  0;
};

/* Split an inclusive integer interval into contiguous blocks with equal square-root
 * increments.  Principal CM representations have density roughly 1/sqrt(|D|), so this
 * balances expected curve orders substantially better than equal-width D intervals. */
scsqrtpart(lo, hi, slot, slots) = {
  if(hi < lo, return([1, 0]));
  my(a = sqrt(lo), b = sqrt(hi+1), first, last);
  first = floor((a + (b-a)*slot/slots)^2);
  last = floor((a + (b-a)*(slot+1)/slots)^2) - 1;
  [max(lo, first), min(hi, last)];
};

/* Split an inclusive interval into equal-width blocks.  This balances the
 * discriminant tests themselves and is preferable for a factor-free first pass. */
scwidthpart(lo, hi, slot, slots) = {
  if(hi < lo, return([1, 0]));
  my(count = hi-lo+1, first, last);
  first = lo + (count*slot) \ slots;
  last = lo + (count*(slot+1)) \ slots - 1;
  [first, last];
};

/* CM-first root search.  Worker slot k of slots scans disjoint contiguous blocks of
 * odd and/or even fundamental discriminants.  Contiguous blocks avoid locking a worker
 * into a residue class with a local obstruction to Cornacchia representations.  Each
 * principal representation 4p=t^2+|D|v^2 supplies the exact candidate orders p+1+-t
 * without SEA.  kind=0 scans both discriminant types, 1 only odd, and 2 only even. */
shortcertcm(p, {slot = 0}, {slots = 1}, {dstart = 3}, {dbound = 100000}, {smoothbits = 40}, {kind = 0}, {reset = 1}, {partition = 0}) = {
  if(!ispseudoprime(p), error("short: p is composite"));
  if(p < 5, error("short: need p >= 5"));
  if(slot < 0 || slot >= slots, error("short: invalid CM worker slot"));
  if(kind < 0 || kind > 2, error("short: invalid CM discriminant kind"));
  if(partition < 0 || partition > 1, error("short: invalid CM partition"));
  my(n = #binary(p), n2 = n^2, L = scbound(p), rt = sqrtint(p),
     d, d0, kmin, kmax, bounds, firstk, lastk, d0min, d0max, first0, last0,
     v, cert, oldrootsmooth = SC_rootsmoothbits);
  if(reset, screset(p));
  SC_rootsmoothbits = smoothbits;
  if(kind != 2,
    kmin = max(0, dstart \ 4);
    kmax = (dbound-3) \ 4;
    bounds = if(partition, scwidthpart(kmin, kmax, slot, slots),
                scsqrtpart(kmin, kmax, slot, slots));
    firstk = bounds[1]; lastk = bounds[2];
    forstep(d = 3 + 4*firstk, 3 + 4*lastk, 4,
      if(SC_maxcurves && SC_curves >= SC_maxcurves, break);
      SC_cmtests++; SC_cmlastD = d;
      if(quaddisc(-d) != -d, next);
      /* A representation 4p=t^2+d*v^2 makes -d a square modulo p.
       * The Kronecker symbol is much cheaper than the modular square root
       * attempted by qfbcornacchia, and rejects about half the remaining
       * fundamental discriminants without losing any representation. */
      if(kronecker(-d, p) == -1, next);
      v = qfbcornacchia(d, 4*p);
      if(!#v, next);
      cert = sccmtrace(p, n2, L, rt, -d, v[1]);
      if(type(cert) == "t_VEC",
        SC_rootsmoothbits = oldrootsmooth;
        return(cert))
    )
  );
  if(kind != 1,
    d0min = max(1, (dstart+3) \ 4);
    d0max = dbound \ 4;
    bounds = if(partition, scwidthpart(d0min, d0max, slot, slots),
                scsqrtpart(d0min, d0max, slot, slots));
    first0 = bounds[1]; last0 = bounds[2];
    for(d0 = first0, last0,
      if(SC_maxcurves && SC_curves >= SC_maxcurves, break);
      SC_cmtests++; SC_cmlastD = 4*d0;
      if(quaddisc(-d0) != -4*d0, next);
      /* Likewise, p=t^2+d0*v^2 implies that -d0 is a square modulo p. */
      if(kronecker(-d0, p) == -1, next);
      v = qfbcornacchia(d0, p);
      if(!#v, next);
      cert = sccmtrace(p, n2, L, rt, -4*d0, 2*v[1]);
      if(type(cert) == "t_VEC",
        SC_rootsmoothbits = oldrootsmooth;
        return(cert))
    )
  );
  SC_rootsmoothbits = oldrootsmooth;
  0;
};

/* Resume a checkpointed CM root order (stored with xden=0).  Factoring is repeated only
 * on the smallest saved residual; the class polynomial is still deferred until the order
 * passes the certificate-window test. */
shortcertcmfromorder(p, D, N, {R0 = 0}, {reset = 1}) = {
  if(!ispseudoprime(p), error("short: p is composite"));
  if(p < 5, error("short: need p >= 5"));
  if(D >= 0 || quaddisc(D) != D, error("short: invalid CM discriminant"));
  my(n = #binary(p), n2 = n^2, L = scbound(p), rt = sqrtint(p),
     C, lev, tail, childstop);
  if(reset, screset(p));
  SC_resumeattempts++;
  C = scscreenorder(p, n2, N, L, rt, D, R0);
  if(C == 0, return(0));
  sclogcmscreen(p, D, N, C);
  lev = scmontgomerylevel(p, D, N, C);
  if(type(lev) != "t_VEC",
    if(lev == 0, sclogcmtombstone(p, D, N, C));
    return(0));
  scloglevel(p, lev);
  if(lev[4] != 1, SC_descents++);
  childstop = if(SC_branchcurves, SC_curves + SC_branchcurves, 0);
  tail = scchain(lev[4], n2, childstop);
  if(type(tail) != "t_VEC", SC_backtracks++; return(0));
  concat([p, lev[1], lev[2], lev[3]], tail);
};

/* Resume a CM root that has already passed factoring and window selection.  The saved
 * (D,N,o,q) tuple avoids repeating the discriminant scan and root-order factorization;
 * reconstruction and the exact-order point test are still performed from scratch. */
shortcertcmfromscreen(p, D, N, o, q, {reset = 1}) = {
  if(!ispseudoprime(p), error("short: p is composite"));
  if(p < 5, error("short: need p >= 5"));
  if(D >= 0 || quaddisc(D) != D, error("short: invalid CM discriminant"));
  my(n = #binary(p), n2 = n^2, L = scbound(p), rt = sqrtint(p),
     t = abs(N-p-1), w, v2, sr, m, fo, C, lev, tail, childstop);
  if(N % 2, error("short: saved CM curve order is odd"));
  w = 4*p - t^2;
  if(w <= 0 || w % (-D), error("short: saved CM order has an invalid trace"));
  v2 = w / (-D);
  if(sqrtint(v2)^2 != v2, error("short: saved CM order has an invalid norm"));
  if(N % o, error("short: saved CM certificate order does not divide the curve order"));
  if(q == 1,
    m = o,
    if(q <= n2 || q >= rt || !ispseudoprime(q) || o % q,
      error("short: invalid saved CM child prime"));
    m = o/q
  );
  if(reset, screset(p));
  SC_resumeattempts++;
  sr = smoothpart(o, n2);
  if(sr[2] != q, error("short: saved CM child does not match certificate order"));
  m = sr[1];
  if(m <= 1 || o <= L || o >= factor(m)[1,1] * L,
    error("short: saved CM order is outside the certificate window"));
  fo = factor(m)[,1];
  if(q != 1, fo = concat(fo, [q]~));
  C = [[o, q, fo]];
  lev = scmontgomerylevel(p, D, N, C);
  if(type(lev) != "t_VEC",
    if(lev == 0, sclogcmtombstone(p, D, N, C));
    return(0));
  scloglevel(p, lev);
  if(lev[4] != 1, SC_descents++);
  childstop = if(SC_branchcurves, SC_curves + SC_branchcurves, 0);
  tail = scchain(lev[4], n2, childstop);
  if(type(tail) != "t_VEC", SC_backtracks++; return(0));
  concat([p, lev[1], lev[2], lev[3]], tail);
};

/* Resume an exact root curve order saved by SC_candidatefile.  If F is supplied, it must
 * be a row or column of known factors of the saved rough residual; otherwise the currently
 * configured bounded factoring portfolio is applied again.  SEA is not recomputed. */
shortcertfromorder(p, A, xden, N, {F = 0}, {R0 = 0}, {reset = 1}) = {
  if(!ispseudoprime(p), error("short: p is composite"));
  if(p < 5, error("short: need p >= 5"));
  my(n = #binary(p), n2 = n^2, L = scbound(p), rt = sqrtint(p),
     d = scnonsquare(p), E, sr, s, dv, lev, tail, childstop);
  if(reset, screset(p));
  SC_resumeattempts++;
  if(xden == 1,
    E = ellinit([0, A, 0, 1, 0], p),
    if(xden != d, error("short: invalid saved twist denominator"));
    E = ellinit([0, (d*A)%p, 0, (d^2)%p, 0], p)
  );
  if(#E == 0, error("short: saved curve is singular"));
  if(type(F) == "t_VEC" || type(F) == "t_COL",
    if(type(F) == "t_VEC", F = F~);
    sr = smoothpart(N, n2); s = sr[1];
    SC_maxsmoothbits = #binary(s);
    dv = divisors(s);
    lev = sctryfactors(p, n2, A, xden, E, N, L, rt, s, dv, F),
    lev = sctryorder(p, n2, A, xden, E, N, L, rt, R0)
  );
  if(lev == 0, return(0));
  scloglevel(p, lev);
  if(lev[4] != 1, SC_descents++);
  childstop = if(SC_branchcurves, SC_curves + SC_branchcurves, 0);
  tail = scchain(lev[4], n2, childstop);
  if(type(tail) != "t_VEC", SC_backtracks++; return(0));
  concat([p, lev[1], lev[2], lev[3]], tail);
};

/* Resume a complete root level that already passed the exact-order point test.  Cheap
 * structural checks reject malformed checkpoints; the final verifier still validates
 * the curve and point before the runner accepts any completed certificate. */
shortcertfromlevel(p, A, x, o, q, {reset = 1}) = {
  if(!ispseudoprime(p), error("short: p is composite"));
  if(p < 5, error("short: need p >= 5"));
  my(n = #binary(p), n2 = n^2, L = scbound(p), sr, m, tail, childstop);
  if(reset, screset(p));
  SC_resumeattempts++;
  if(A < 0 || A >= p || gcd(A^2-4, p) != 1,
    error("short: invalid saved curve parameter"));
  if(x < 0 || x >= p || o < 1, error("short: invalid saved root level"));
  sr = smoothpart(o, n2); m = sr[1];
  if(sr[2] != q, error("short: saved child does not match root order"));
  if(m <= 1 || o <= L || o >= factor(m)[1,1] * L,
    error("short: saved root order is outside the certificate window"));
  if(q != 1 && (q <= n2 || q^2 >= p || !ispseudoprime(q)),
    error("short: invalid saved child prime"));
  if(q == 1, return([p, A, x, o]));
  SC_descents++;
  childstop = if(SC_branchcurves, SC_curves + SC_branchcurves, 0);
  tail = scchain(q, n2, childstop);
  if(type(tail) != "t_VEC", SC_backtracks++; return(0));
  concat([p, A, x, o], tail);
};

printshort(p) = {                                        \\ "p A0 x0 o0 A1 x1 o1 ..."
  my(c = shortcert(p));
  for(i = 1, #c, printf("%d%s", c[i], if(i < #c, " ", "\n")));
};
