set datafile separator ';'
set xrange [-5:5]
set yrange [-5:5]
#set zrange [0:1]
#splot 'neccomp.csv' using 1:2:(20*log10($3)) w p ti "FRICO Ex", 'neccomp.csv' using 7:8:(20*log10($9)) w p ti "NEC Ex"
#splot 'neccomp.csv' using 1:2:(20*log10($4)) w p ti "FRICO Ez", 'neccomp.csv' using 7:8:(20*log10($10)) w p ti "NEC Ez"
#set zrange [0:20]
splot 'neccomp.csv' using 1:2:(20*log10(($3)/($9))) w p ti "difference Ex"
#splot 'neccomp.csv' using 1:2:(20*log10(($4)/($10))) w p ti "difference Ez"
