BEGIN{
  FS = "[ \t]+"
  OFS = "\t"
}
FNR==NR{
  sid = $sc
  s   = tolower($xc)
  if      (s ~ /^(male|m|1)$/)   g = 1
  else if (s ~ /^(female|f|2)$/) g = 2
  else                           g = -9
  sexmap[sid] = g
  next
}
{
  # Use FNR==1 (first line of the *second* file) to detect/preserve header
  if (FNR==1 && ($1=="FID" || $1=="#FID")) {
    print $0, "SEX"
    next
  }
  iid = $2
  g = (iid in sexmap) ? sexmap[iid] : -9
  print $0, g
}
