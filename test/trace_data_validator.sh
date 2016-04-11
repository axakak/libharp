#!/bin/bash

TDT="../bin/traceDataTester"
#INPUTTRACEDATA="traceDataFake.txt"
INPUTTRACEDATA="stc_1_sample/ak_2014113132133.yaml"

$TDT $INPUTTRACEDATA

#TODO: add validation for export file

# sed -i '' -e 's/data/events/g' *.txt     # replace data tag with event tag
# sed -i '' -e '/^[[:space:]]*$/d' *.txt   # delete blank lines
# sed -i '' -e 's/^location: $/&""/' *.txt # add empty quotes after location (note: $ is not working)
# sed -i '' -E 's/^(patient-id:|location:) /&""/' *.txt

#sed -i '' -e 's/data/events/g' -e '/^[[:space:]]*$/d' -E 's/^(patient-id:|location:) /&""/' *.txt

# sed -i '' -e '/workspace:/ i\ # insert classification tag befor workspace tag
# classification:\
# \  group: 0\
# '

#cat trace_list.txt | xargs grep -l "type: width" | xargs sed -i '' -e 's/\(.*group: \)\([[:digit:]]\)$/\11/'

#grep -lr "complexity" * | xargs grep -l "level: 1" >> comp_1_list.txt

# change file ext
# for f in *.txt; do mv $f "${f%.*}.yaml"; done
