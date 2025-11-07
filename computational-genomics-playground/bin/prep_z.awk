# prep_z.awk
# Usage:
#   awk -f prep_z.awk -v BIM=region.pruned.bim -v PIN=prune.in -v OUT=out.tsv -F'[ \t]+' zfile.tsv
#
# - Parses sumstats (tabs OR spaces)
# - Reads BIM to build CHR:BP map
# - For each line in prune list, uses only FIRST TOKEN (handles corrupted lines)
# - Prefer exact SNP id; else fallback to CHR:BP via BIM
# - Writes rows in the same order as prune list

function normchr(x,   y){ y=x; gsub(/^chr/i, "", y); return y }
function poskey(c,b){ return normchr(c) ":" b }
function split_snp_to_pos(s,   a){ n=split(s,a,":"); return (n>=2)? poskey(a[1],a[2]) : "" }

BEGIN{
  OFS="\t"
  # Read BIM → id2pos & pos2id (tab or space delimited)
  while ((getline line < BIM) > 0){
    n = split(line, b, /[ \t]+/)
    if (n >= 4){
      chr=b[1]; id=b[2]; bp=b[4]
      pk = poskey(chr,bp)
      if (id!="") id2pos[id]=pk
      if (pk!="" && !(pk in pos2id)) pos2id[pk]=id
    }
  }
  close(BIM)
}

# First file is the zfile: header + data by SNP id and CHR:BP
FNR==1 {
  snp_col=0
  for (i=1;i<=NF;i++) if ($i=="SNP") snp_col=i
  if (snp_col==0){ print "ERROR: SNP column not found in zfile header" > "/dev/stderr"; exit 1 }
  header=$0
  next
}
{
  s=$snp_col
  if (s!=""){
    line_by_id[s]=$0
    pk=split_snp_to_pos(s)
    if (pk!="") line_by_pos[pk]=$0
  }
  next
}

END{
  print header > OUT
  kept=missing=0; misslist=""

  # stream prune list; use only first token per line
  while ((getline p < PIN) > 0){
    if (p=="") continue
    split(p, t, /[ \t]+/)
    key=t[1]
    if (key=="") continue

    if (key in line_by_id){
      print line_by_id[key] >> OUT
      kept++
    } else {
      pk = (key in id2pos) ? id2pos[key] : split_snp_to_pos(key)
      if (pk!="" && (pk in line_by_pos)){
        print line_by_pos[pk] >> OUT
        kept++
      } else {
        missing++
        if (missing<=5) misslist = (misslist=="" ? key : misslist "," key)
      }
    }
  }
  close(PIN); close(OUT)
  printf("[PREP_Z] kept=%d missing=%d\n", kept, missing) > "/dev/stderr"
  if (missing>0) printf("[PREP_Z] examples missing: %s\n", misslist) > "/dev/stderr"
  if (kept==0) print "[PREP_Z] WARNING: produced only header (no matches)." > "/dev/stderr"
}
