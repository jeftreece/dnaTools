#!/bin/bash

cd $REDUX_ENV
cd /Users/jazdrv/_prj/dnatools.zak/_env/redux2

DB="variant"
TB="$1"

echo ""
sqlite3 $DB.db <<!
.schema $TB
!
echo ""
