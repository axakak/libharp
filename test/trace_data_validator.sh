#!/bin/bash

TDT="../bin/traceDataTester"
#INPUTTRACEDATA="traceDataFake.txt"
INPUTTRACEDATA="stc_1_sample/ak_2014113132133.yaml"

$TDT $INPUTTRACEDATA

#TODO: add validation for export file

#sed -i '' -e 's/data/events/g' *.txt
#sed -i '' -e '/^[[:space:]]*$/d' *.txt
#sed -i '' -e 's/\(^location: ""\).*$/\1/' *
#sed -i '' -e 's/data/events/g' -e '/^[[:space:]]*$/d' -e 's/^patient-id: $/&""/' -e 's/^location: $/&""/' *.txt
#cat stc_2_width_2_list.txt | xargs grep -l "type: width" | xargs sed -i '' -e 's/\(.*group: \)\([[:digit:]]\)$/\11/'

#grep -lr "complexity" * | xargs grep -l "level: 1" >> comp_1_list.txt
