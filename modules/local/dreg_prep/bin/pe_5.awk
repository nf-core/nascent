BEGIN {
	OFS = "\t"
}

($10 == "-") {
	print $1, $6 - 1, $6, $7, $8, $10
}

($10 == "+") {
	print $1, $5, $5 + 1, $7, $8, $10
}
