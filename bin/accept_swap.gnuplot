set multiplot
set size 0.5, 0.9
set origin 0,0

ylow = 0.005
yhi = 1
set key left
set log
set title "mc^3 swap accept prob. run1"
plot [][ylow:yhi] "./x" \
	using ($2/$3) t"T1:T2", \
	"" using ($4/$5) t"T2:T3", \
	"" using ($6/$7) t"T3:T4", \
	"" using ($8/$9) t"T1:T3", \
	"" using ($10/$11) t"T2:T4", \
	"" using ($12/$13) t"T1:T4"

set origin 0.45, 0
set title "mc^3 swap accept prob. run2"
plot [][ylow:yhi] "./x" \
        using ($14/$15) t"T1:T2", \
        "" using ($16/$17) t"T2:T3", \
        "" using ($18/$19) t"T3:T4", \
        "" using ($20/$21) t"T1:T3", \
        "" using ($22/$23) t"T2:T4", \
        "" using ($24/$25) t"T1:T4"
set nomultiplot

pause -1
