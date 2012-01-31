set log 
x_low = 500
y_low = 1e-3

set multiplot
set size 0.5, 1
set origin 0,0
 plot [x_low:*][1:300] "./MB.converge" using 1:6 t"ESS tree length", \
	"" using 1:8 t"ESS alpha", \
	"" using 1:10 t"ESS pInv"

set size 0.5, 1
set origin 0.45,0
plot [x_low:*][y_low:*] \
	"" using 1:($7-1) t"b/w treelength", \
	"" using 1:($9-1) t"b/w alpha", \
	"" using 1:($11-1) t"b/w pInv", \
	"" using 1:5 t"avg splits stddev", \
	"" using 1:12 t"KSD LogL", \
	"" using 1:13 t"KSD TL", \
	"" using 1:14 t"KSD alpha", \
	"" using 1:15 t"KSD pinv", 1000/x

unset multiplot
pause -1
