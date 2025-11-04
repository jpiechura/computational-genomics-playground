# bin/parse_clumped.awk
# Robust parser for PLINK .clumped (space-padded, leading spaces, variable column order)
BEGIN { OFS = "\t" }

# skip empty lines
/^[[:space:]]*$/ { next }

# helper: trim-left + split on any whitespace
function split_ws(line, arr) {
  sub(/^[[:space:]]+/, "", line)
  return split(line, arr, /[[:space:]]+/)
}

# header: build index by column name, then print tidy header
header==0 {
  line = $0
  n = split_ws(line, f)
  for (i=1; i<=n; i++) idx[f[i]] = i
  print "CHR","SNP","BP","P","SP2"
  header = 1
  next
}

# data rows: trim/split and select by names; require CHR to look like 1..22/X/Y/MT
header==1 {
  line = $0
  n = split_ws(line, f)
  if (f[1] ~ /^[0-9XYMT]+$/) {
    # Handle SP2 possibly being missing; default to last field if not present
    sp2_col = ("SP2" in idx ? idx["SP2"] : n)
    print f[1], f[idx["SNP"]], f[idx["BP"]], f[idx["P"]], f[sp2_col]
  }
}
