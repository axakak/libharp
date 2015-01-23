#!/bin/bash

TDT="../bin/traceDataTester"
#INPUTTRACEDATA="traceDataFake.txt"
INPUTTRACEDATA="sample_traces/stc_1/20130828161520.txt"

$TDT $INPUTTRACEDATA

#TODO: add validation for export file

#sed -i '' -e 's/data/events/g' *.txt
