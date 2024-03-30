#!/bin/awk -f
BEGIN {
	OFS = "\t"
}

{
	print $1, $2, $3, -1 * $4
}

