#!/bin/bash

. ~/Code/github/evsel/bshm/evsel.bshm

evflow |\
 evrangeyear "${1}" "${2}" |\
 evrangemag "${3}" "${4}" |\
 evrangedepth "${5}" "${6}" |\
 evrangelat "${7}" "${8}" |\
 evrangelon "${9}" "${10}" |\
 evprepare lon lat depth mag > ${11}
 