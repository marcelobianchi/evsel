#!/bin/bash

##
# This is a wrapper code used by evgraph. You should not use it
# alone. There is support neither documentation for this !
##

[ -z "$EVSELBASE" ] && EVSELBASE=/usr/local/share/evsel/
[ ! -d $EVSELBASE/ ] && echo "Cannot find BSHM file @ '$EVSELBASE', please define EVSELBASE correctly."
. $EVSELBASE/evsel.bshm

evflow |\
 evrangeyear "${1}" "${2}" |\
 evrangemag "${3}" "${4}" |\
 evrangedepth "${5}" "${6}" |\
 evrangelat "${7}" "${8}" |\
 evrangelon "${9}" "${10}" |\
 evprepare lon lat depth mag > ${11}
 