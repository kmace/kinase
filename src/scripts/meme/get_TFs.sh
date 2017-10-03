cat $1 | grep "\"alt\":" | awk -F"\"" '{print $4}'
