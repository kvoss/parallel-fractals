#!/bin/tcsh

#foreach opt ( "-m" "-o" "")
foreach opt ( "-s4" "-s8")
foreach np ( 1 2 4 8 )
	foreach ds ( 1 2 4 8 )
		@ dss = $ds * 1000

		@ ttime = 0
		foreach i ( `seq 3` )
			set t = `./mandelbrot_set -n $np -r $dss $opt | grep 'Elapsed' | awk '{print $3*1000}'`
			@ ttime += $t
		end
		@ ttime /= 3

		echo "numproc: $np, resolution: $ds [k], options: $opt, time: $ttime [ms]"
	end
end
end
