#!/bin/bash

cd $REDUX_PATH
cd /Users/jazdrv/_prj/dnatools/_env/redux2

DB="variant"
TB="$1"
LIM="$2"

if [ "$TB" == "" ]; then
    echo "Need a tbl as first arg. Aborting."
    exit
fi
if [ "$LIM" == "" ]; then
    LIM=5
fi

echo ""
sqlite3 $DB.db <<!
.header on
.mode column
select count(*) from $TB;
select * from $TB limit $LIM;
!
echo ""
