#!/bin/awk -f
BEGIN {
	OFS = "\t"
}

{
	print $1, $2, $3, $4 * 1000 * 1000 / readCount / 1
}
