
run:
	Rscript -e "shiny::runApp('../shinyGMM')"

lib:
	R CMD SHLIB "src/read_bin.c"

packages:
	Rscript mypkg.R
