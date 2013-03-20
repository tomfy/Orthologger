set terminal x11 enhanced size 900,300
set multiplot
set size 0.33, 0.99
set origin 0,0

ylow = 0.005
yhi = 1
set key left
set log
set title "mc^3 swap accept prob. run1"
plot [][ylow:yhi] "./x" \
	using 1:($2/$3) t"T1:T2", \
	"" using 1:($4/$5) t"T2:T3", \
	"" using 1:($6/$7) t"T1:T4" 
set origin 0.33, 0
set title "mc^3 swap accept prob. run2"
plot [][ylow:yhi] "" \
     using 1:($8/$9) t"T1:T2", \
	"" using 1:($10/$11) t"T2:T3", \
	"" using 1:($12/$13) t"T1:T3"

set origin 0.67, 0
set title "mc^3 swap accept prob. run3"
plot [][ylow:yhi] "./x" \
        using 1:($14/$15) t"T1:T2", \
        "" using 1:($16/$17) t"T2:T3", \
        "" using 1:($18/$19) t"T1:T3"
#, \
      #  "" using 1:($20/$21) t"T1:T3", \
#        "" using 1:($22/$23) t"T2:T4", \
#        "" using 1:($24/$25) t"T1:T4"
set nomultiplot

pause -1
