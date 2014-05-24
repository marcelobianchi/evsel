#!/bin/bash 

status=0
for tool in awk cat wget rm sort wc mv sed mkdir
do
	which $tool > /dev/null 2>&1
	[ $? -ne 0 ] && echo "Failed $tool" && status=1
done
exit $status

for tool in disaz
do
	which $tool > /dev/null 2>&1
	[ $? -ne 0 ] && echo "$tool was not found, some functions may fail"
done

exit $status
