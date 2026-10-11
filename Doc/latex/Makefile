# Extracted source archive: no Python package or network required.
BOOKS = main theory application api verification
.PHONY: all
all:
	@set -e; for book in $(BOOKS); do \
	  pdflatex -no-shell-escape -halt-on-error -interaction=nonstopmode $$book.tex; \
	  bibtex $$book; makeindex $$book.idx; \
	  pdflatex -no-shell-escape -halt-on-error -interaction=nonstopmode $$book.tex; \
	  pdflatex -no-shell-escape -halt-on-error -interaction=nonstopmode $$book.tex; \
	done
