set multiplot
set size 0.5, 0.9
set origin 0,0


plot [][0:*] "./x" \
	using ($2/$3) t"Tdiff 1 level", \
	"" using ($4/$5) t"Tdiff 1 level", \
	"" using ($6/$7) t"Tdiff 1 level", \
	"" using ($8/$9) t"Tdiff 2 levels", \
	"" using ($10/$11) t"Tdiff 2 levels", \
	"" using ($12/$13) t"Tdiff 3 levels"

set origin 0.45, 0
plot [][0:*] "./x" \
        using ($14/$15) t"Tdiff 1 level", \
        "" using ($16/$17) t"Tdiff 1 level", \
        "" using ($18/$19) t"Tdiff 1 level", \
        "" using ($20/$21) t"Tdiff 2 levels", \
        "" using ($22/$23) t"Tdiff 2 levels", \
        "" using ($24/$25) t"Tdiff 3 levels"
set nomultiplot

pause -1
