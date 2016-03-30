ssh patator -t "cd $1; echo \"knitr::knit('$2')\" | R --vanilla"
