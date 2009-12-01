#!/bin/sh

TEMPDIR=`mktemp -d`

## Download pycircuit git repo
git clone git://github.com/henjo/pycircuit.git $TEMPDIR/pycircuit

## Generate docs
cd $TEMPDIR/pycircuit
ls
. ./setup.sh
cd doc
if make html | tee $TEMPDIR/logfile
then
    ## Update web documentation on success
    cp -fR $TEMPDIR/pycircuit/doc/build/html/* /nfs/data/www/docs.pycircuit.org
else
   ## Email logfile if build failed
    true
fi

rm -fr $TEMPDIR