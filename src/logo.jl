using Luxor

Drawing(625,625,"logo.svg")

origin()
background("white")
setfont("Lithos Pro Black", 100)
N,F,C = juliacircles()
settext("F",F+Point(0,12),halign="center",valign="center")
settext("N",N+Point(0,12),halign="center",valign="center")
settext("C",C+Point(0,12),halign="center",valign="center")

finish()
preview()