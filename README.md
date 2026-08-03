# ShortPrimalityProofs
For the purpose of this repository, a **short ECPP** is a sequence of integers $(p_0,A_0,x_0,m_0p_1,A_1,x_1,\ldots,m_{k-1}p_k,A_k,x_k,m_kp_{k+1})$ in which
- we have $k \ge 0$, $p_0 \ge 5$, and put $n:=\lceil \log_2 p_0\rceil$,
- the $p_i$ are odd integers that satisfy $n^2 < p_{i+1} < \sqrt{p_i}$ for $0\le i < k$ and $p_{k+1}=1$,
- the $m_i$ are $n^2$-smooth integers satisfying $L_i < m_ip_{i+1} < r_iL_i$, where $L_i = q_i+1+\lfloor 2\sqrt{q_i}\rfloor$ with $q_i=\lfloor\sqrt{p_i}\rfloor$ and $r_i$ is the least prime divisor of $m_i$,
- each $A_i$ is a nonnegative integer less than $p_i$ with $\gcd(A_i^2-4,p_i)=1$,
- each $x_i$ is a nonnegative integer less than $p_i$,

such that for each $0\le i \le k$ there exist integers $B_i,y_i\in [0,p_i-1]$ with $\gcd(B_i,p_i)=1$ for which $(x_i,y_i)$ is a point of order $m_ip_{i+1}$ on the [Montgomery curve](https://en.wikipedia.org/wiki/Montgomery_curve) $B_iy^2 = x^3 + A_ix^2 +x$ modulo every prime divisor of $p_i$.

For integers $p_0 > 3$, a short ECPP $(p_0,...)$ exists if and only if $p_0$ is prime.  The size of a short ECPP is linear in $n=\log p_0$ and can be verified in quasi-quadratic time (in contrast to the quadratic size and quasi-cubic verification time for a conventional ECPP).

This repository contains the following resources:
- vsmallECPP.py is a Python program that verifies a short ECPP in quasi-quadratic time.
- short8all.txt contains the 201,072 short ECPPs with $p_0\le 2^8$.
- short.gp is a GP script that uses either exact CM orders or SEA on random curves and their quadratic twists, with bounded backtracking from unlucky descents, to search for short ECPPs.
- parallel_short.py runs independent short.gp searches with a portfolio of complete factoring, partial factoring, and SEA early-abort strategies, and verifies the first certificate found.
- certs.csv is a list of short ECPPs for the primes listed in the table below.

For example, to search for a certificate with four workers, write the first verified result to `cert.txt`, and append durable run metadata to `search-runs.jsonl`, run
```
python3 parallel_short.py $(echo 'nextprime(10^150)' | gp -q) -j 4 -o cert.txt --manifest search-runs.jsonl
```

The manifest records the prime, start and finish times, number of parallel worker attempts,
every worker's seed and search configuration, wall and aggregate child CPU time, the winning
worker, and the last curve/SEA counters reported by each worker.  Explicit nonconsecutive seeds
can be supplied with `--seeds`.  With ten workers, the default portfolio uses eight complete-factor
workers and two partial-factor workers.  Experimental PARI `ellsea` torsion filters can be assigned
per worker with `--sea-torsions`; the default is unfiltered point counting because rejecting orders
with unwanted small factors does not necessarily improve certificate-search throughput.
Counters remain cumulative when a worker moves from a resumed order or CM presearch into its
random-SEA fallback, so the manifest accounts for all work performed by that attempt.
For large searches, `--deep-factor-bits` and `--deep-factor-seconds` can reserve a longer
factorization budget for an order whose rough cofactor has already been reduced below a chosen
bit length, without spending that budget on every order.  This adaptive pass can be used even
when the worker's regular factor budget is zero.
Setting a worker's `--factor-seconds` entry to zero disables its final PARI `factorint` pass while
retaining saved external GMP-ECM work, which is useful when short timed PARI calls repeatedly
discard their partial progress.
`--prefactor-seconds` adds a short saved partial pass before the full factorization; flag 8 is the
default for this stage because it omits PARI's expensive final ECM step while retaining factors
already found.  Every returned candidate is still tested and any accepted descent is recursively
proved, so this does not weaken certificate verification.
The GP search enables PARI's `factor_add_primes` facility for every timed `factorint` call.  After
an alarm expires, a cheap bounded factorization recovers prime factors that PARI discovered before
the timeout, reduces and checkpoints the residual, and reports the event as `factor_recoveries`.
Recovery uses `factor(R,2)`: since `R` is already $n^2$-rough, this consults PARI's add-prime
table without repeating trial division through a large bound.
This applies to random SEA orders, exact CM orders, resumed roots, and every recursive level.
`--factor-rounds` bounds automatic continuation (three rounds by default); a further round is run
only when the preceding pass found factors and strictly reduced the residual.
All saved-factor stages stop as soon as `smooth * min(residual, sqrt(p))` cannot reach the
certificate lower bound, since no prime factor left in that residual can then produce a descent.
Optional `--pm1-bounds` and `--pp1-bounds` prepasses use GMP's P-1 and P+1 methods before ECM;
any factors they find are retained and the reduced cofactor is passed to the later stages.
`--ecm-rounds` likewise permits another GMP-ECM batch only after the preceding batch found a
factor and reduced the residual; a no-factor batch is never repeated automatically.
When `msieve` and GNU `timeout` are installed, `--msieve-bits` and `--msieve-seconds` add a
one-thread, wall-time-bounded msieve stage for residuals below a selected size; every factor it
emits before completion or timeout is retained.
The exact $n^2$-smooth part of every curve order is extracted with repeated gcds against a cached
product of all primes through $n^2$.  This replaces a separate division by every small prime for
each order while retaining every prime power exactly; at the 210-digit target the cache is under
100 KiB and is built once per worker.  Curve family 5 implements the optimized $X_1(27)$ model,
retaining Montgomery-compatible curves with full rational 2-torsion and a rational point of order
4, so the selected-side group order is divisible by 216 before point counting.
The $X_1(27)$ finite-field model and map are adapted from
[`IslayResearch/OneShotSEA`](https://github.com/IslayResearch/OneShotSEA), which implements and
authenticates Andrew Sutherland's published optimized $X_1(27)$ formula; license attribution is
retained in `THIRD_PARTY_NOTICES.md`.

An optional CM-first root search avoids point counting for its candidate orders.  For a negative
fundamental discriminant $D$, Cornacchia's algorithm can cheaply detect a representation
$4p=t^2+|D|v^2$, which gives the two exact candidate orders $p+1\pm t$.  The search trial-divides
these orders first and spends the configured P-1, P+1, ECM, prefactor, msieve, and bounded
`factorint` stages only when the smooth part reaches `--cm-smooth-bits` (40 by default).  Saved
factors are checkpointed exactly as they are for SEA-derived orders.  A usable order is converted
back to the same Montgomery certificate format through its Hilbert class polynomial; all recursive
levels and final verification are unchanged.  Since this Montgomery model has the rational
2-torsion point $(0,0)$, odd CM group orders are impossible and are discarded before any factoring.
`--cm-start A --cm-bound B` partitions the inclusive discriminant
range into contiguous, square-root-weighted worker blocks before each worker falls back to
its ordinary random-SEA search.  This weighting balances the expected number of Cornacchia
representations, whose density decreases roughly as $1/\sqrt{|D|}$, much better than equal-width
blocks;
`--cm-kind odd`, `even`, or `both` selects the fundamental-discriminant families.  Using
non-overlapping ranges makes consecutive batches reproducible and avoids duplicated work.
For a factor-free enumeration pass, `--cm-partition width` instead gives each worker the same
number of discriminants.  This balances wall time when testing discriminants dominates; the
default `sqrt` partition remains useful when inline factoring cost follows the expected number
of represented orders.
Pass `--cm-only` to make this a finite range scan: if every worker exhausts its assigned block,
the runner exits successfully and writes an `outcome: "exhausted"` finish record instead of
starting the random-SEA fallback.

Long searches can preserve promising root orders instead of discarding them after a bounded
factorization.  For example, `--candidate-bits 48 --candidate-dir candidates/180` writes each
order whose $n^2$-smooth part has at least 48 bits to a per-worker append-only file.  A record is
`p A xden N s R`: the prime, curve parameter, twist denominator, exact curve order, smooth part,
and remaining composite residual.  CM orders use `xden=0` and store $|D|$ in the `A` field, so
they use the same checkpoint files without constructing a curve first.  Repeated records for one
order checkpoint smaller residuals found by saved factoring stages.  Once an exact-order point
test succeeds, the same file also receives a five-field `p A x o q` record for the complete root
level.  This preserves the most expensive progress even if its recursive child subtree later
backtracks.  A CM order that passes factoring and window selection is saved before class-polynomial
reconstruction as a seven-field `p |D| N o q 0 0` record, so an interrupted reconstruction does
not lose the viable order.  Passing the file or directory back with `--resume-candidates` ranks
complete root levels first, then screened CM roots, followed by rough saved orders ranked by the
estimated size of the complementary cofactor (then by smoothness), assigns them to workers, and
resumes without repeating SEA or root factoring before falling back to fresh search.
`--resume-top N` limits a portfolio to the strongest N checkpoints.  The GP function
`shortcertfromorder(p,A,xden,N)` provides the same resume path directly; an optional row or column
of known residual factors can be supplied as its fifth argument.  Passing `0` there and a saved
residual as the sixth argument resumes directly from that checkpoint after verifying that it
divides the order's original rough part.  `shortcertfromlevel(p,A,x,o,q)` resumes a complete
five-field root-level checkpoint; the final verifier still validates the entire resulting proof.
`shortcertcmfromscreen(p,D,N,o,q)` resumes a seven-field screened CM root and repeats only model
reconstruction, the exact-order point test, and the recursive child proof.
`--resume-per-worker N` lets each worker consume a round-robin queue of N globally ranked
checkpoints in one process.  Counters and factor tables remain live across the queue, avoiding
repeated startup and allowing a two-pass CM workflow: first enumerate exact orders with all
factoring limits set to zero, then apply the expensive portfolio only to the strongest saved
orders.
For example, the first command below performs a balanced finite enumeration without factoring,
and the second spends the full portfolio on the top 200 resulting checkpoints, twenty per worker:
```
python3 parallel_short.py "$p" -j 10 --cm-bound 100000000 --cm-screen-only \
  --candidate-bits 40 --candidate-dir candidates/210 --manifest search-runs-210.jsonl
python3 parallel_short.py "$p" -j 10 --resume-candidates candidates/210 \
  --resume-top 200 --resume-per-worker 20 --resume-only --candidate-bits 40 \
  --candidate-dir candidates/210 --manifest search-runs-210.jsonl -o cert210.txt
```
The checkpoint loader rechecks that five- and seven-field orders have exactly one admissible
$n^2$-rough child and that their complementary factor is $n^2$-smooth before ranking them.
After a class polynomial is constructed successfully but supplies no matching Montgomery
level, the runner marks each attempted seven-field screen as exhausted so that it is skipped
on later runs.  Once every known screen for an order is exhausted, its rough checkpoint is
also suppressed.  Reconstruction errors do not write tombstones and remain retryable.
Use `--resume-only` for a finite ranked-candidate pass.  The runner then records an exhausted
outcome after all assigned candidates fail, rather than silently continuing into CM or random SEA.
Fully factored but unusable orders are written back with residual `1` as tombstones; the loader also
discards any order whose entire saved residual is below the exact minimum admissible child prime,
so an earlier, larger checkpoint for an already exhausted order is not replayed.
`--resume-offset K --resume-top N` selects a disjoint ranked slice for staged portfolios.

For large inputs, install PARI's optional
[`seadata` package](https://pari.math.u-bordeaux.fr/packages.html).  It supplies modular-polynomial
tables used by `ellcard`/`ellsea` and covers the inputs in this challenge through $10^{200}$;
PARI uses these tables automatically for every SEA call, including recursive levels and all
parallel workers.  `parallel_short.py` reports their availability at startup, records it in the
run manifest, and warns when they are unavailable.  The search remains correct without them, but
falls back to substantially slower modular-polynomial computations.

**Challenge**

Below is a list of short ECPPs for the least prime $p>10^c$ for $c=10,20,\ldots,200$.
The entries through $c=100$ were found by <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a>
on a single core.  The entries from $c=110$ through $c=200$ were found by <a href="parallel_short.py">parallel_short.py</a>,
which runs single-threaded short.gp workers in parallel; the listed CPU time is for the winning
worker unless the entry is explicitly labeled as a checkpointed, combined total.

<details>
<summary>$p=10^{10}+19$,&nbsp; via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (&lt;1 CPU second, 2 levels).</summary>

```
10000000019 9322349340 1921958667 116108 7217 235 607
```
</details>
<details>
<summary>$p=10^{20}+39$,&nbsp; via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~1 CPU seconds, 2 levels).</summary>

```
100000000000000000039 89951393186720294033 26135327929659638076 11876954936 509599 419481 3373
```
</details>
<details>
<summary>$p=10^{30}+57$,&nbsp; via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~1 CPU seconds, 2 levels).</summary>

```
1000000000000000000000000000057 106991342299430347297585638871 188675017969977395328028624483 9586166844055967 14077590431530 22429269288900 4934867
```
</details>
<details>
<summary>$p=10^{40}+121$,&nbsp; via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~2 CPU seconds, 3 levels).</summary>

```
10000000000000000000000000000000000000121 4902559719197972567355860705693483269960 1108182596968252581482904098615171436619 136850847522421485837 4858820576325 5128710645172 5627651 10966 13232 267
```
</details>
<details>
<summary>$p=10^{50}+151$,&nbsp; via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~3 CPU seconds, 3 levels).</summary>

```
100000000000000000000000000000000000000000000000151 4168493225324236537663200519121316886619997624556 45590400574778487393338352639345746341371438526095 13601869893526828282016090 333961040995 170482403318 935683 40410 98277 1193
```
</details>
<details>
<summary>$p=10^{60}+7$,&nbsp; via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~3 CPU seconds, 3 levels).</summary>

```
1000000000000000000000000000000000000000000000000000000000007 582594560733647942864167554683101932559817757617443321218616 485107000290978194104100671498974843211121245518759155383025 1591141197950235962428006613308 13059006973177308093777996339 13059586929974475337549368522 1179275706066947 10409174197 23529769559 248836
```
</details>
<details>
<summary>$p=10^{70}+33$,&nbsp; via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~45 CPU seconds, 4 levels).</summary>

```
10000000000000000000000000000000000000000000000000000000000000000000033 9193638238367761016751675951961267177328031093306760076142420392370661 5154446894383162465732893340634268054473714420328923353004710402276202 177664059059812136202711786749417049 8649705917060546335363234070525383 42605680596938635666326245539183459 303885057908347530 9987604800409038 5847492647777907 1449221531 3881343 7598118 3432
```
</details>
<details>
<summary>$p=10^{80}+129$,&nbsp; via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~180 CPU seconds, 4 levels).</summary>

```
100000000000000000000000000000000000000000000000000000000000000000000000000000129 21738970501572142787334250995459457313706731047873015853200788703257646728746387 41991271598373420878627810537757671546773541549398474815901013810480424900203390 10475053568008371237180368561737331009607 6251483180616232563867118006 2989838890071853567090370520 109497485294994 64373892028 116329932035 967431233 82783 99666 1889
```
</details>
<details>
<summary>$p=10^{90}+289$,&nbsp; via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~350 CPU seconds, 3 levels).</summary>

```
1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000289 292390063605290292583797761883066555291812969778434353799388301927489872556173779693570186 201936337656195409865372057596506482368540069004779260915505630735608561247615266092698548 1495231769958646021642577303700692519195356608 1370055486155935701230818919713 631674234039720505250594205971 2824893046457853 65743 243188 50291
```
</details>
<details>
<summary>$p=10^{100}+267$,&nbsp; via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~800 CPU seconds, 3 levels).</summary>

```
10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000267 9922776908302697864551893522633781539036740348871476332653810465626471659307742084010027213820226625 2836186773449440978191705565856816079355043173337102496065404967247933485306358946184015108636857285 159990589984239373593745787852818347617199852606512 989772819515561949937579106033692378969575949199 160869327784807806952128639954291424038666088760 2979396535159070660937591 221054722104 68038539734 625635
```
</details>
<details>
<summary>$p=10^{110}+7$,&nbsp; via <a href="parallel_short.py">parallel_short.py</a> (~850 CPU seconds, 4 levels).</summary>

```
100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000007 6463512687231274453314405725127389227686510737885438111434665718277396637427547966320233024608351230230861450 94164843423036637004338287352695725396243376738721874335250825798680411803544771226508486391278261756148716303 11887982306821438660772844495374643156319942615511664684 202249972180221490910665359620639420738815558212027 150718343273118010945877719422581588244803222920063 69271420438922496826047729 1783250918118183 7920143168078823 192942401 8303705 3827745 37339
```
</details>
<details>
<summary>$p=10^{120}+79$,&nbsp; via <a href="parallel_short.py">parallel_short.py</a> (~1400 CPU seconds, 4 levels).</summary>

```
1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000079 132063756553348037607726246459351520545782379718387545073787885595088443428910224053262081547377506070464921665690244443 766759587882597019696391814500801893207393762520787020231259750746628420452950574472606593873044215447309872888211933399 2032407269668174614737384431875199029029184284197878747552629 1199197362586136210273591905540461080971671185565224891285 1235746649075380226211220902889073034006519549905970133438 103755726490907726500461806304 3092901622776628753789127 800991492200408803139269 9138053287529 5714221 2420272 45971
```
</details>
<details>
<summary>$p=10^{130}+1113$,&nbsp; via <a href="parallel_short.py">parallel_short.py</a> (~55 CPU seconds, 3 levels).</summary>

```
10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001113 8838054562218699739469769827167347178366084538631527476350638177688407070489981665391546772154550105442112460156690474042240671153 9070442852984080644726137002845969709600170053521757279077458871926798567371260193902248343382297986025014308252012522028148277940 143023612271167604986843854406903897909541795703955557853382262527 4279311677016680608246200050783064414358575781642362807565916 2514922654280625968927317951358588049284674438860504597506761 4852063280663853797388579475379 8756592747341738788031594 5272134364285153528036270 3037506078068
```
</details>
<details>
<summary>$p=10^{140}+13$,&nbsp; via <a href="parallel_short.py">parallel_short.py</a> (~800 CPU seconds, 4 levels).</summary>

```
100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000013 58641989056675720896672212411665741567276120290417855224286073593918069783779639337562924504400725658476104789609335447371195336334882769838 49477426077173743212373028968808725584389507235326871104842053964179679773766627681933870798077570405431407409265070868732921824050415203356 13102807158286493606695640041882330489449163784474568128908900189019724 9194143353218612803141831381890939699022032316343544916837584131946 9446062608280649498405523563298491468876161342077046283467239660178 3772692773649837804777542519049033 104014555298455253904175742770 56869381954496753251959470195 815915526665845 1388267 2031435 164771
```
</details>
<details>
<summary>$p=10^{150}+67$,&nbsp; via <a href="parallel_short.py">parallel_short.py</a> (~1500 CPU seconds, 5 levels).</summary>

```
1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000067 318242351116513938119517967616301475671486714041557198261085687482277145984397246332451109797710555221454953956380848003465751737440305039178134859211 720477903566384043117699259849645626560674563858596916881593754258176690470197473175283690734334846930259881683458635858051341569534187208064885697331 3274688382283250370121362465461297838165277819878561512256860826179463643421 419710093855186023983988396848350352203085299959848344457841346392366791 2738868097477462328198525091918736337577217055081977682817874517862422856 7790490239009640360339423052989945581 165044629407820523237731842524831306 65747092211780530324214107412979080 536903698434235914 19738996086605 20356708900543 4294067863 719611 318088 127969
```
</details>
<details>
<summary>$p=10^{160}+303$,&nbsp; via <a href="parallel_short.py">parallel_short.py</a> (~1200 CPU seconds, 5 levels).</summary>

```
10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000303 1459610408796319836784355481875424878795873921791031879135324980121367679324525174627177955281717548199131578517243540832793275941249177793495795195459013374976 4758106707451115267653405186080802099032088051348450733120020318362579889791887278021953210724205680736180811567428047248223479702468810707973314657809819169814 1019270628659222703537437424625804232803263645948203567236028599943201735592927793 448568730017408308873597708912961564403152659118218528354017946397296103365610 43320014997290486975014588988392937147131394524080329101272755177864009581337 700117564993251844165020303999846974095 718795565085543471731276397799125154 10968362732951735055925268568007215413 6744096990366654554 2535320783744295006 618160122705634570 4589520461 900519 465044 1420
```
</details>
<details>
<summary>$p=10^{170}+63$,&nbsp; via <a href="parallel_short.py">parallel_short.py</a> (~830 CPU seconds, 4 levels).</summary>

```
100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000063 97851301468710322096701415255025735445020416732173420584513156348596460580522182207203840964259659298440287710679350631819479721427397620850835974604321152741859027353616 55385901170942206534243384784986217946125130313590212796246065260973505826112767900263571516970685154718325092034578736234215098064253126247743096031225951704979818277124 14313191454343359675134982037110351112302714967446571694914676917005824035072641255903 48290115752735916780161082915022062724389380253590340158786626525365801795262948 41400976291934482046519894194588855766774648118129282571744903347019212194151736 12014709593774384763720483892757011363852 1354285071734168727488863597074859438 1475043512004542153372361051227787964 1457813217634777254 3131982 1160249 101809
```
</details>
<details>
<summary>$p=10^{180}+313$,&nbsp; via <a href="parallel_short.py">parallel_short.py</a> (~850 CPU seconds, 4 levels).</summary>

```
1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000313 676314487666624426613264577413941673599493662940655898562057070611245460551288155535091786832776359738580085216599663216771952967980765637893778558559819437935500208875865390249244 437588306890387722589892662578409413226705638943666450657740601038982508421822919939552462866951504536517925265322210790219141179699027744923317252688881621374183616123925848009969 1166336987321735003266566202037904971880547636635894093000838667118924529020328958360209285 181230328020199888788167879924035658953885968492820660612531541546092140880882001291413109 180740652885697327946185328655571807826790152811195448489571030474428954342453352279735131 594606096126489693565799487615702466460829756 347840253512180047137361416902967552057 985111460480449331112129520819167682825 87106255631465576389 282050657021841 118686513408715 58134903
```
</details>
<details>
<summary>$p=10^{190}+253$,&nbsp; via <a href="parallel_short.py">parallel_short.py</a> (~1400 CPU seconds, 4 levels).</summary>

```
10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000253 9032059704630233730649942650536047579693026386223857885968948358385423848143230060965273666340984028744349337999958813471343224132869169824548586931328844850141135248586236633667573758081519 8672697054356251821949201159081436227181229023012836043139512670476451145681061017199275263310290576408282862769512107642454036464908582668082645931949598601721146074781572407560952323413337 119560819481213251567972071174401671086589625916629105075129472983232261018104163860364995031733 210153835907840010764977747557600768073435083199840813626221127092220892700368741636592 124131114470027993484417042859632405674168231719628843745417560746617765125842629893875 69380085617604829723575538604964918163436573 755176061281954099938074618632834941680 529594975764018078624919429248545500806 5089982853083610347687 12937680819 21605446100 309513
```
</details>
<details>
<summary>$p=10^{200}+357$,&nbsp; via <a href="parallel_short.py">parallel_short.py</a> (~4400 combined CPU seconds across the saved root and resumed child chain, 4 levels).</summary>

```
100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000357 57834291042150292495995563609652096282941917786582520751588842815371340027216658579432999024002709915305656128581620025390448845210225124177561529171174256999204857514501328789133046559706504654019024 11014890427879667372499868798052112552439497346445545565585187250033094195651188657919359377356221017943370375110895658588552650470341345570751912995867508645166912655969863486486372574381190051284073 450153423952056675364935857477723043101127369120523333039334166374620323711184334771772552020494112201 403197499065195869777358579896162701416242962142916541281892395320601106094138547554226079864635718 2289902501950831138158126198253297541276141821724838373701307085709523139401377355655500975010297949 87313221154928125278237123659305289785138559754195 669697371797517697784437372418971920470508070459 338761438897941287847471585101516158723344790380 2887645374971118626289461 4732035594544 19377378581309 9920677
```
</details>
