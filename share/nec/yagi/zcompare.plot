plot 'zdata.txt' using 1:2 w lp t 'R', \
     'zdata.txt' using 1:3 w lp t 'X', \
     'yagi.gplot' using ($1)*1e6:2 w lp t 'Rnec', \
     'yagi.gplot' using ($1)*1e6:3 w lp t 'Xnec'
