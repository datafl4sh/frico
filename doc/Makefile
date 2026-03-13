all: frico.html frico.pdf

frico.html: ../man/frico.1
	cat ../man/frico.1 | groff -mandoc -Thtml > frico.html

frico.pdf: ../man/frico.1
	cat ../man/frico.1 | groff -mandoc -Tpdf > frico.pdf


