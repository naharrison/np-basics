# https://wiki.sagemath.org/art?action=show&redirect=pics

path_op = dict(rgbcolor='red', thickness=1)
fill_op = dict(rgbcolor='blue', alpha=0.3)
options = dict(pathoptions=path_op, filloptions=fill_op, endarrow=False, startpoint=False)
G = [words.fibonacci_tile(i).plot(**options) for i in range(8)]
a = animate(G)
a.show(delay=150)
