using Luxor

Drawing(625,625,"logo.svg")

origin()
background("transparent")
Luxor.scale(1.7)
setfont("Lithos Pro Black", 180)
N,F,C = juliacircles()
settext("F",F+Point(0,12),halign="center",valign="center")
settext("N",N+Point(0,12),halign="center",valign="center")
settext("C",C+Point(0,12),halign="center",valign="center")

finish()

Drawing(625,625,"logo-dark.svg")

origin()
background("transparent")
Luxor.scale(1.7)
setfont("Lithos Pro Black", 180)
N,F,C = juliacircles()
setcolor("white")
settext("F",F+Point(0,12),halign="center",valign="center")
settext("N",N+Point(0,12),halign="center",valign="center")
settext("C",C+Point(0,12),halign="center",valign="center")

finish()
preview()