## Reference values for the regression tests in test-regression-*.R.
##
## These were computed with excursions 2.5.11.9000 before the rewrite of the
## Takahashi recursion (Qinv), the integration routine (shapeInt) and several
## R functions for speed, so the tests check that the rewrite did not change
## any results. The problems are defined in helper-regression.R, and each
## entry is a summary of one computation, as described in the tests.
##
## The entries for excursions() and continuous() (ex>, ex<, ex=, ex!=, ex.qc
## and cont.*) were recomputed after a fix of the reordering of the nodes: the
## two nodes with the largest marginal probabilities were in the same
## constraint set for CAMD, and are now integrated in the order of their
## probabilities like the other nodes. This changed the results by less than
## 1e-4. The entry sc was recomputed after a.marginal and b.marginal were
## swapped back in simconf(), and its bounds a and b are unchanged.
##
## Monte Carlo results depend on the random number streams, so only change
## these values if a change in the results is intended.

REF <- list(
  gi.natural = c(P = 0.43802827674496803, E = 0.0062263572981268497, Pv.sum = 70.723810093526296, 
    Pv.wsum = 4108.8829918547399, Ev.sum = 0.486807974232213, Ev.wsum = 17.875968224964701
    ),
  gi.mu = c(P = 0.43802827674496803, E = 0.0062263572981268497, Pv.sum = 70.723810093526296, 
    Pv.wsum = 4108.8829918547399, Ev.sum = 0.486807974232213, Ev.wsum = 17.875968224964701
    ),
  gi.maxsize = c(P = 0, E = 0, Pv.sum = 37.377282863121799, Pv.wsum = 3048.1019476187298, 
    Ev.sum = 0.077890573585978506, Ev.wsum = 5.2587733148340901),
  gi.lim = c(P = 0, E = 0, Pv.sum = 67.069195788134095, Pv.wsum = 4092.18615990851, 
    Ev.sum = 0.436015830969954, Ev.wsum = 17.645768818312899),
  gi.ind = c(P = 0.579631042560944, E = 0.0070160633187534198, Pv.sum = 73.084841924902605, 
    Pv.wsum = 4170.1183457412599, Ev.sum = 0.47316747249877, Ev.wsum = 16.215123923595399
    ),
  gi.limits = c(P = 0.60360076049814504, E = 0.0055968807081320602, Pv.sum = 62.659395042088597, 
    Pv.wsum = 3268.8327449017602, Ev.sum = 0.53602493725977896, Ev.wsum = 25.974738804912899
    ),
  gi.sparsity = c(P = 0.443409382719888, E = 0.0052560166933792096, Pv.sum = 60.791713179127797, 
    Pv.wsum = 3019.0696559913199, Ev.sum = 0.44255971578473602, Ev.wsum = 21.859130145489399
    ),
  gi.chol = c(P = 0.43802827674496803, E = 0.0062263572981268497, Pv.sum = 70.723810093526296, 
    Pv.wsum = 4108.8829918547399, Ev.sum = 0.486807974232212, Ev.wsum = 17.875968224964701
    ),
  gi.seed1 = c(P = 0.44625256785614598, E = 0.0063346536740123, Pv.sum = 71.474892763994006, 
    Pv.wsum = 4133.9214197880901, Ev.sum = 0.49571987359605701, Ev.wsum = 18.212086503858799
    ),
  gi.onesided = c(P = 0.67564497997029904, E = 0.0065732487052806198, Pv.sum = 83.497380867412502, 
    Pv.wsum = 4504.62038568761, Ev.sum = 0.45242880801879398, Ev.wsum = 16.300687967357199
    ),
  rand = c(0.0010094978404174399, 0.59500378387998498, 0.35783453761357398, 0.22234082670111499, 
    0.46682759725957601),
  `ex>` = c(F.sum = -80.323347050437505, F.wsum = -4328.4839794572699, E.sum = 10, 
    E.wsum = 368, M.sum = -80, M.wsum = -4314, vars.sum = 22.791760883106999,
    vars.wsum = 1150.9839245969099),
  `ex<` = c(F.sum = -38.5700368517276, F.wsum = -2117.5168367456599, E.sum = 31, 
    E.wsum = 1483, M.sum = -69, M.wsum = -3567, vars.sum = 22.791760883106999,
    vars.wsum = 1150.9839245969099),
  `ex=` = c(F.sum = -64.302057571094593, F.wsum = -3495.1664481149701, E.sum = 35, 
    E.wsum = 1523, M.sum = -59, M.wsum = -3312, vars.sum = 22.791760883106999,
    vars.wsum = 1150.9839245969099),
  `ex!=` = c(F.sum = -30.6979424289053, F.wsum = -2035.8335518850299, E.sum = 35, 
    E.wsum = 1523, M.sum = -59, M.wsum = -3312, vars.sum = 22.791760883106999,
    vars.wsum = 1150.9839245969099),
  ex.ind = c(F.sum = -88.202702207730198, F.wsum = -4265.9261689837404, E.sum = 6, 
    E.wsum = 399),
  ex.qc = c(F.sum = -86.270327636919205, F.wsum = -4601.8524778537803, E.sum = 7, 
    E.wsum = 229),
  cm = c(F.sum = -46.1910278836502, F.wsum = -2474.3591403600999, E.sum = 5, E.wsum = 92, 
    M.sum = -86, M.wsum = -4733, P0 = NA, P1 = 0.73721821938533505, P2 = 0.033960032812370998,
    P0b = 0.68597594343530599, P1b = 0.96734444531691299, P2b = 0.82771169123704202
    ),
  sc = c(a.sum = -42.582893096383003, a.wsum = -780.69186561747904, b.sum = 42.687217565436498, 
    b.wsum = 796.80518162618102, am.sum = -24.469111313470201, am.wsum = -445.58690263359301
    ),
  exmc = c(F.sum = 17.492000000000001, F.wsum = 733.85599999999999, E.sum = 9, E.wsum = 290, 
    M.sum = -82, M.wsum = -4470),
  exmc.ind = c(F.sum = 64.804000000000002, F.wsum = 3215.9859999999999, E.sum = 23, 
    E.wsum = 1210, M.sum = -74, M.wsum = -3655),
  cmmc = c(F.sum = 15.022, F.wsum = 563.05200000000002, E.sum = 5, E.wsum = 92, 
    M.sum = -86, M.wsum = -4733, P1 = 0.76000000000000001, P2 = 0.029999999999999999
    ),
  scmc = c(a.sum = -140.66999465706701, a.wsum = -7076.1137917753404, b.sum = 142.86747168458601, 
    b.wsum = 7363.9220171408097, am.sum = -76.731621146534494, am.wsum = -3813.0446442364801,
    bm.sum = 76.387380483202605, bm.wsum = 3906.7845736077102),
  scmc.ind = c(a.sum = -92.689890258060103, a.wsum = -2810.9468494457801, b.sum = 93.529234153426998, 
    b.wsum = 3012.8174934969202),
  tc.plus = c(nloc = 618, loc1 = 285.87994960819299, loc2 = 316.06554468742701, loc.w1 = 90962.199762039701, 
    loc.w2 = 107607.569034228, nidx = 620, idx1 = 67272585, idx2 = 67637456,
    grp1 = 0, grp2 = 60, grp3 = 5, grp4 = 141, grp5 = 10, grp6 = 220, grp7 = 5,
    grp8 = 179, grp.w = 1088564),
  tc.minus = c(nloc = 445, loc1 = 208.85461889154399, loc2 = 223.330344012171, loc.w1 = 48219.287696675303, 
    loc.w2 = 55278.857245715502, nidx = 447, idx1 = 25490772, idx2 = 25700045,
    grp1 = 0, grp2 = 60, grp3 = 8, grp4 = 187, grp5 = 10, grp6 = 180, grp7 = 2,
    grp.w = 455662),
  tc.onlevel = c(nloc = 279, loc1 = 117.926216320084, loc2 = 153.48078663925401, loc.w1 = 15169.9752760421, 
    loc.w2 = 26031.8663492054, nidx = 284, idx1 = 6794547, idx2 = 6865835,
    grp1 = 0, grp2 = 50, grp3 = 5, grp4 = 92, grp5 = 8, grp6 = 125, grp7 = 4,
    grp.w = 191569),
  seg1 = c(n = 3, len1 = 5, len2 = 4, len3 = 3, seq1 = 35, seq2 = 67, seq3 = 62, 
    seg1 = 30, seg2 = 37, seg3 = 25, grp1 = 17, grp2 = 6, grp3 = 9),
  seg2 = c(n = 4, len1 = 3, len2 = 4, len3 = 3, len4 = 3, seq1 = 10, seq2 = 63, 
    seq3 = 58, seq4 = 14, seg1 = 4, seg2 = 35, seg3 = 26, seg4 = 11, grp1 = 3,
    grp2 = 6, grp3 = 9, grp4 = 6),
  seg3 = c(n = 4, len1 = 78, len2 = 61, len3 = 101, len4 = 19, seq1 = 677586, seq2 = 451817, 
    seq3 = 1159300, seq4 = 50454, seg1 = 649232, seg2 = 456082, seg3 = 1257014,
    seg4 = 37590, grp1 = 11631, grp2 = 3660, grp3 = 20200, grp4 = 684),
  cont.step = c(F.sum = 22.093875436561198, F.wsum = 3109.5259125662201, M = 112, Msum = 56.3333333333333, 
    Mw = 3230.5555555555602),
  cont.linear = c(F.sum = 35.778768513192702, F.wsum = 5718.3966037660302, M = 192, Msum = 93.630651539414401, 
    Mw = 9858.3864816036803),
  cont.log = c(F.sum = 22.331549129995601, F.wsum = 3156.4153417207899, M = 112, Msum = 56.3333333333333, 
    Mw = 3230.5555555555602)
)

## Integration with very narrow intervals for some nodes, computed as above
REF$gi.narrow <- c(
  P = 4.3289074339864702e-46, E = 6.1125633705659002e-47,
  Pv.sum = 5.0002558394319498, Pv.wsum = 490.02314952654899,
  Ev.sum = 3.9434179075516802e-06, Ev.wsum = 0.00035835680349186798
)

## Reference values for simconf.mixture. These were computed after the
## mixture quantiles were changed to be computed to full precision instead of
## with the default tolerance of uniroot, which changed the bands by about
## 3e-5, and after the fixes of the sequential integration method.
REF.MIXTURE <- list(
  int = c(
    a.sum = -93.237583627468695, a.wsum = -3008.3179910620402,
    b.sum = 89.057498587154399, b.wsum = 2916.2721809132199,
    am.sum = -47.322809394417298, am.wsum = -1516.08782848786,
    bm.sum = 50.427094578328997, bm.wsum = 1660.7840506263899
  ),
  samp = c(
    a.sum = -92.936089124539706, a.wsum = -2998.5194197168398,
    b.sum = 88.807544094138095, b.wsum = 2908.1486598901802,
    am.sum = -47.322809394417298, am.wsum = -1516.08782848786,
    bm.sum = 50.427094578328997, bm.wsum = 1660.7840506263899
  )
)
