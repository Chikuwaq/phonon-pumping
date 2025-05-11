CoThicknessNm = 36.1
PtThicknessNm = 21.32

set term pdfcairo enhanced color font "Helvetica,18" size 5in,3in
set output sprintf("FMR_CoThicknessNm%.1f_PtThicknessNm%.1f.pdf", CoThicknessNm, PtThicknessNm)
set lmargin screen 0.20
set rmargin screen 0.80

set pm3d map
set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000') # define color gradation. Relative to cbrange

set title sprintf("FMR power absorption: Co/Pt %.1f/%.1f nm", CoThicknessNm, PtThicknessNm)

set xlabel 'Out-of-plane magnetic field [T]'
set ylabel 'Frequency [GHz]'

sp 'spectrum.dat' notitle w pm3d

unset pm3d

unset term
set output