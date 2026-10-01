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
## The entries for the contours (tc.*, seg3 and cont.*) summarise the
## geometry in a way that does not depend on the numbering of the mesh, which
## differs between platforms. See the fingerprint functions in
## helper-regression.R.
##
## All Monte Carlo entries (gi.*, ex*, cm, sc, cont.* and REF.MIXTURE$int)
## were recomputed after the sequential integration was changed to draw the
## samples in chunks with one random stream for each chunk, which makes the
## results independent of the number of threads. The new values agree with
## the old ones up to the Monte Carlo error.
##
## Monte Carlo results depend on the random number streams, so only change
## these values if a change in the results is intended.

REF <- list(
  gi.natural = c(P = 0.44776099271419689, E = 0.0062719345474585428, Pv.sum = 71.506785685765237, 
    Pv.wsum = 4135.3129882664061, Ev.sum = 0.48888450839098063, Ev.wsum = 17.927012956232527
    ),
  gi.mu = c(P = 0.44776099271419689, E = 0.0062719345474585428, Pv.sum = 71.506785685765237, 
    Pv.wsum = 4135.3129882664061, Ev.sum = 0.48888450839098063, Ev.wsum = 17.927012956232527
    ),
  gi.maxsize = c(P = 0, E = 0, Pv.sum = 37.43605047177261, Pv.wsum = 3051.8234752627427, 
    Ev.sum = 0.077888842204821895, Ev.wsum = 5.2514892442366872),
  gi.lim = c(P = 0, E = 0, Pv.sum = 67.768939396015185, Pv.wsum = 4118.2341648008214, 
    Ev.sum = 0.43761454890361928, Ev.wsum = 17.694520270799824),
  gi.ind = c(P = 0.59260871328542819, E = 0.0070437441050739315, Pv.sum = 73.940329559064139, 
    Pv.wsum = 4198.8228718601622, Ev.sum = 0.47546566611648466, Ev.wsum = 16.303083214203678
    ),
  gi.limits = c(P = 0.58628882692552753, E = 0.0056502098046811826, Pv.sum = 61.067601127250263, 
    Pv.wsum = 3194.5987523515405, Ev.sum = 0.54074240844227539, Ev.wsum = 26.189426224228384
    ),
  gi.sparsity = c(P = 0.44481845577478685, E = 0.0053657342753381469, Pv.sum = 60.967501608765374, 
    Pv.wsum = 3028.3432251868508, Ev.sum = 0.44529238395403564, Ev.wsum = 22.092685173905679
    ),
  gi.chol = c(P = 0.44776099271419689, E = 0.006271934547458541, Pv.sum = 71.506785685765252, 
    Pv.wsum = 4135.312988266407, Ev.sum = 0.48888450839098052, Ev.wsum = 17.927012956232527
    ),
  gi.seed1 = c(P = 0.43262065361385416, E = 0.0062217822827513817, Pv.sum = 70.381976597733242, 
    Pv.wsum = 4096.9268530249674, Ev.sum = 0.48863616562960849, Ev.wsum = 17.934024211420859
    ),
  gi.onesided = c(P = 0.69132021152364742, E = 0.0066383399730184454, Pv.sum = 84.447431894714285, 
    Pv.wsum = 4535.5240991787969, Ev.sum = 0.45140471549463074, Ev.wsum = 16.222527349282554
    ),
  rand = c(0.0010094978404174399, 0.59500378387998498, 0.35783453761357398, 0.22234082670111499, 
    0.46682759725957601),
  `ex>` = c(F.sum = -80.328686996673213, F.wsum = -4328.611876162091, E.sum = 10, 
    E.wsum = 368, M.sum = -80, M.wsum = -4314, vars.sum = 22.791760883107049, 
    vars.wsum = 1150.9839245969058),
  `ex<` = c(F.sum = -38.552667737449987, F.wsum = -2116.6119244476868, E.sum = 31, 
    E.wsum = 1483, M.sum = -69, M.wsum = -3567, vars.sum = 22.791760883107049, 
    vars.wsum = 1150.9839245969058),
  `ex=` = c(F.sum = -64.333708777212522, F.wsum = -3496.5537605967183, E.sum = 35, 
    E.wsum = 1523, M.sum = -59, M.wsum = -3312, vars.sum = 22.791760883107049, 
    vars.wsum = 1150.9839245969058),
  `ex!=` = c(F.sum = -30.666291222787486, F.wsum = -2034.4462394032817, E.sum = 35, 
    E.wsum = 1523, M.sum = -59, M.wsum = -3312, vars.sum = 22.791760883107049, 
    vars.wsum = 1150.9839245969058),
  ex.ind = c(F.sum = -88.200204806302466, F.wsum = -4265.7454700854478, E.sum = 6, 
    E.wsum = 399),
  ex.qc = c(F.sum = -86.271379990978758, F.wsum = -4601.8392899596192, E.sum = 7, 
    E.wsum = 229),
  cm = c(F.sum = -47.477648938301606, F.wsum = -2505.8886674770715, E.sum = 5, 
    E.wsum = 92, M.sum = -86, M.wsum = -4733, P0 = NA, P1 = 0.73659049253652431, 
    P2 = 0.031576332089226992, P0b = 0.68597594343530632, P1b = 0.96734444531691299, 
    P2b = 0.82771169123704247),
  sc = c(a.sum = -42.488462809393617, a.wsum = -778.94490530817666, b.sum = 42.592787278447197, 
    b.wsum = 795.05822131687819, am.sum = -24.469111313470201, am.wsum = -445.58690263359301
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
  tc.plus = c(nloc = 618, loc1 = 285.87994960819299, loc2 = 316.06554468742701, loc3 = 127741.71468111299, 
    loc4 = 97034.486880877099, nidx = 620, edges1 = 287.12842375069101, edges2 = 317.67325376595301,
    edges3 = 285.86807463555903, edges4 = 316.70222924456101, edges5 = 3536,
    edges6 = 128751.038857169, edges7 = 97730.456333036404, edges8 = 127745.748065765,
    edges9 = 98453.700240614897, edges10 = 1126038, grp1 = 0, grp2 = 60, grp3 = 5,
    grp4 = 141, grp5 = 10, grp6 = 220, grp7 = 5, grp8 = 179),
  tc.minus = c(nloc = 445, loc1 = 208.85461889154399, loc2 = 223.330344012171, loc3 = 66573.475183506394, 
    loc4 = 48516.503648151302, nidx = 447, edges1 = 210.42353706421699, edges2 = 224.33217896181301,
    edges3 = 209.765405775787, edges4 = 223.06789675353099, edges5 = 2036,
    edges6 = 67328.588634891305, edges7 = 48843.698519199497, edges8 = 66872.713441708998,
    edges9 = 49333.729430458101, edges10 = 463032, grp1 = 0, grp2 = 60, grp3 = 8,
    grp4 = 187, grp5 = 10, grp6 = 180, grp7 = 2),
  tc.onlevel = c(nloc = 279, loc1 = 117.926216320084, loc2 = 153.48078663925401, loc3 = 24570.590994507998, 
    loc4 = 20504.401158415702, nidx = 284, edges1 = 120.426216320084, edges2 = 156.62220799549101,
    edges3 = 118.61524315084201, edges4 = 155.33936528301601, edges5 = 1301,
    edges6 = 25598.629343135799, edges7 = 21346.893734054702, edges8 = 24974.0046063874,
    edges9 = 21646.984781091502, edges10 = 190801, grp1 = 0, grp2 = 50, grp3 = 5,
    grp4 = 92, grp5 = 8, grp6 = 125, grp7 = 4),
  seg1 = c(n = 3, len1 = 5, len2 = 4, len3 = 3, seq1 = 35, seq2 = 67, seq3 = 62, 
    seg1 = 30, seg2 = 37, seg3 = 25, grp1 = 17, grp2 = 6, grp3 = 9),
  seg2 = c(n = 4, len1 = 3, len2 = 4, len3 = 3, len4 = 3, seq1 = 10, seq2 = 63, 
    seq3 = 58, seq4 = 14, seg1 = 4, seg2 = 35, seg3 = 26, seg4 = 11, grp1 = 3,
    grp2 = 6, grp3 = 9, grp4 = 6),
  seg3 = c(n = 4, seqs1 = 259, seqs2 = 4, seqs3 = 892, seqs4 = 106.284177414244, 
    seqs5 = 128.78149111978399, seqs6 = 779, seqs7 = 10, seqs8 = 2812, seqs9 = 303.20430430792101,
    seqs10 = 387.93083215784401),
  cont.step = c(F1 = 180.50000000000009, F2 = 180.49999999999972, F.F = 22.077847295583041, 
    F4 = 44102.166666666664, F5 = 33272.166666666657, F.F = 4301.7064029426283, 
    nrings = 4, rings1 = 1, rings2 = 52, rings3 = 1.0370370370370374, rings4 = 27.888888888888886, 
    rings5 = 24.888888888888886, rings6 = 4, rings7 = 148, rings8 = 3.0864197530864206, 
    rings9 = 81.444444444444429, rings10 = 74.222222222222214),
  cont.linear = c(F1 = 180.50000000000009, F2 = 180.49999999999972, F.F = 35.758968660231488, 
    F4 = 44102.166666666664, F5 = 33272.166666666657, F.F = 7065.4909235135101, 
    nrings = 4, rings1 = 1, rings2 = 92, rings3 = 1.0451357462083386, rings4 = 49.952524299511325, 
    rings5 = 40.342803692826806, rings6 = 4, rings7 = 246, rings8 = 3.1030634394485199, 
    rings9 = 142.78318650972423, rings10 = 119.38732403736188),
  cont.log = c(F1 = 180.50000000000009, F2 = 180.49999999999972, F.F = 22.320119086677501, 
    F4 = 44102.166666666664, F5 = 33272.166666666657, F.F = 4354.6414941134808, 
    nrings = 4, rings1 = 1, rings2 = 52, rings3 = 1.0370370370370374, rings4 = 27.888888888888886, 
    rings5 = 24.888888888888886, rings6 = 4, rings7 = 148, rings8 = 3.0864197530864206, 
    rings9 = 81.444444444444429, rings10 = 74.222222222222214)
)

## Integration with very narrow intervals for some nodes, computed as above
REF$gi.narrow <- c(P = 7.4842837455803139e-46, E = 1.4563977453597287e-46, Pv.sum = 5.0002454768548699, 
  Pv.wsum = 490.02221030671251, Ev.sum = 4.0239367193973504e-06, Ev.wsum = 0.00036625531820737978
  )

## Reference values for simconf.mixture. These were computed after the
## mixture quantiles were changed to be computed to full precision instead of
## with the default tolerance of uniroot, which changed the bands by about
## 3e-5, and after the fixes of the sequential integration method.
REF.MIXTURE <- list(
  int = c(a.sum = -93.222261540777865, a.wsum = -3007.8200232445824, b.sum = 89.044794991900332, 
    b.wsum = 2915.8593140674589, am.sum = -47.32280939441727, am.wsum = -1516.0878284878629, 
    bm.sum = 50.427094578328976, bm.wsum = 1660.7840506263899),
  samp = c(
    a.sum = -92.936089124539706, a.wsum = -2998.5194197168398,
    b.sum = 88.807544094138095, b.wsum = 2908.1486598901802,
    am.sum = -47.322809394417298, am.wsum = -1516.08782848786,
    bm.sum = 50.427094578328997, bm.wsum = 1660.7840506263899
  )
)
