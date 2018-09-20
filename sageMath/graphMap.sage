# http://www.sagemath.org/tour-graphics.html

stnc = 'I am a cool multiedge graph with loops'
g = DiGraph({}, loops=True, multiedges=True)
for a,b in [(stnc[i], stnc[i+1]) for i in xrange(len(stnc)-1)]:
   g.add_edge(a, b, b)
g.plot(color_by_label=True, edge_style='solid').show(figsize=(8,8))
