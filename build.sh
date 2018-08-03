#!/bin/sh
#
# Simple script for building book. Add flag --final to run bibtex,
# makeindex etc and build the final version of the book.

FILE=$1

CH_DIR=chapter_build

# rebuild the entire book
if [ "x$2" = "x--final" ]; then
    pdflatex -shell-escape $FILE;
    bibtex $FILE;
    makeindex $FILE;
    pdflatex -shell-escape $FILE;
    pdflatex -shell-escape $FILE;

# build just a specific chapter :
elif [ "x$2" = "x--chapter" ]; then
		cp $FILE.bbl "$CH_DIR/$3.bbl";
    pdflatex -shell-escape \
              -jobname=$CH_DIR/$3 "\includeonly{tex/$3}\input{$FILE}";
    pdflatex -shell-escape \
             -jobname=$CH_DIR/$3 "\includeonly{tex/$3}\input{$FILE}";
    # remove the parts we had to make to compile the chapter :
		rm $CH_DIR/$3.aux $CH_DIR/$3.bbl $CH_DIR/$3.idx $CH_DIR/$3.log;
		rm $CH_DIR/$3.out $CH_DIR/$3.run.xml $CH_DIR/$3-blx.bib;

# or just recompile :
else
    pdflatex -shell-escape $FILE;
fi
