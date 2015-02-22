#!/bin/bash

TDT="../bin/traceDataTester"
#INPUTTRACEDATA="traceDataFake.txt"
INPUTTRACEDATA="sample_traces/stc_1/20130828161520.txt"

$TDT $INPUTTRACEDATA

#TODO: add validation for export file

#sed -i '' -e 's/data/events/g' *.txt
#sed -i '' -e '/^[[:space:]]*$/d' *.txt
#sed -i '' -e 's/data/events/g' -e '/^[[:space:]]*$/d' -e 's/^patient-id: /&""/' -e 's/^location: /&""/' *.txt
