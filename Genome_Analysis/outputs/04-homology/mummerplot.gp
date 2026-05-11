set terminal png tiny size 800,800
set output "mummerplot.png"
set ytics ( \
 "tig00000005" 1.0, \
 "tig00000006" 40013.0, \
 "tig00000008" 55027.0, \
 "tig00000003" 79862.0, \
 "tig00000009" 108734.0, \
 "tig00000004" 125216.0, \
 "tig00000001" 139949.0, \
 "tig00000007" 2915080.0, \
 "tig00000002" 2931201.0, \
 "" 3147208 \
)
set size 1,1
set grid
unset key
set border 10
set tics scale 0
set xlabel "CP111853.1"
set ylabel "QRY"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
set xrange [1:2803703]
set yrange [1:3147208]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "mummerplot.fplot" title "FWD" w lp ls 1, \
 "mummerplot.rplot" title "REV" w lp ls 2
