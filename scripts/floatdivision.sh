#!/bin/bash

function add_division_import() {
    if ! grep "from __future__ import division" $1 > /dev/null; then
	if tr "\n" " " < $1 | sed  "s/\"\"\".*\"\"\"//g" | sed s/\#.*$//g | sed  "s/'.*'//g" | sed  "s/\".*\"//g" | grep [^/]/[^/] >/dev/null
	then
	    echo Adding division import to $1
	    tmpfile=`mktemp`
	    (echo from __future__ import division;  cat $1) > $tmpfile
	    cat $tmpfile >$1
	    rm -f $tmpfile
	fi
    fi
}

export -f add_division_import

find $1 -name "*.py" -exec bash -c "add_division_import {}" \;
