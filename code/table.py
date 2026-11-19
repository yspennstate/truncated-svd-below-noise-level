import csv, sys, math, statistics as st
from collections import OrderedDict
for path in sys.argv[1:]:
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    g = OrderedDict()
    for r in rows:
        g.setdefault((r['c'], r['T_name']), []).append(r)
    print(path.split('/')[-1], 'rows', len(rows), 'agree', sum(r['agree'] in ('1', 'True', 'true') for r in rows))
    for (c, tn), rs in g.items():
        L = [float(r['loss']) for r in rs]; p = float(rs[0]['pred_op'])
        gap = max(abs(x - p) / p for x in L)
        v2 = [float(r['meas_v2']) for r in rs if r['meas_v2'] not in ('', 'nan')]
        v2 = [x for x in v2 if not math.isnan(x)]
        pv2 = rs[0]['pred_v2']
        se = lambda v: st.stdev(v) / math.sqrt(len(v)) if len(v) > 1 else float('nan')
        s2 = ('%.4f +- %.4f / %s' % (st.mean(v2), se(v2), pv2[:6])) if v2 else '--'
        print('c=%s %-14s T=%s n=%d kept=%s loss %.4f +- %.4f pred %.4f gap %.3f | v2 %s' % (
            c, tn, rs[0]['T'][:6], len(rs), rs[0]['kept'], st.mean(L), se(L), p, gap, s2))
