for dir in `find . -type d`; do [ "$dir" != "." ] && mkdir -p $BASILISK/ibm/$dir ;done
for file in `find . -type f`; do ln -s $IBM_SRC/$file $BASILISK/ibm/$file; done
