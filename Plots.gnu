#input should be formatted as:
#gnuplot -e "systemwidth=...; plots=...; highlight=...; tfinal=..." Plots.gnu
#where plots should be smaller than 10000

unset key
set term png

set title 'Position correlation function of argon (MD simulated with LJ potential)'
set xlabel 'x in units of sigma'
set ylabel 'normalised <x(0)x(r)>'
set output './Data/correlation.png'
plot './Data/correlation.md' with lines

set xrange[0 to tfinal]

set title 'Energy fluctuations in the argon system (MD simulated with LJ potential)'
set xlabel 'time in units of 10^-14s'
set ylabel 'energy in units of epsilon'
set output './Data/energy.png'
plot './Data/energy.md' using 1:2 with lines, './Data/energy.md' using 1:3 with lines, './Data/energy.md' using 1:4 with lines

set title 'Velocity correlation function of argon (MD simulated with LJ potential)'
set xlabel 'time in units of 10^-14s'
set ylabel 'normalised <v(0)v(t)>'
set output './Data/velocityautocorrelation.png'
plot './Data/velocitycorr.md' using 1:2 with lines

unset xrange
unset xlabel
unset ylabel
unset zlabel

set border 4095
set ticslevel 0
unset tics

yview = 45

set xrange[0 to systemwidth]
set yrange[0 to systemwidth]
set zrange[0 to systemwidth]

do for [i=0:plots]{
	set view 60,yview
	
	if(i <= 9){
		outfile = sprintf('./Data/Plot000%i.png',i)
	}
	if(i >= 10 && i <= 99){
		outfile = sprintf('./Data/Plot00%i.png',i)
	}
	if(i >= 100 && i <= 999){	
		outfile = sprintf('./Data/Plot0%i.png',i)
	}
	if(i >= 1000){
		outfile = sprintf('./Data/Plot%i.png',i)	
	}
	set output outfile
	
	set title 'Motion of argon atoms (MD simulated with LJ potential)'	

	splot './Data/coordinates.md' index i every::0::highlight-1 using 3:4:5 pointtype 7 lt rgb "#FF0000", './Data/coordinates.md' index i every::highlight::highlight using 3:4:5 pointtype 7 lt rgb "#00FF00", './Data/coordinates.md' index i every::highlight+1 using 3:4:5 pointtype 7 lt rgb "#FF0000" 
	yview = yview - 0.2
	if(yview < 0){
		yview = yview + 360
	} 	
}
