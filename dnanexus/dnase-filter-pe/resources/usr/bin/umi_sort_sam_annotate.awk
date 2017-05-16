BEGIN{
  FS="\t";
  OFS="\t";
}
{
  if ( $9 > 0 ) {
    if( and($2, 64) ) { 
      strand = "+"; 
    } else { 
      strand = "-";
    }

    match($0, /XD:Z:([AGCT]+)/, umi_match);

    print $1, $5, $3, $4 - 1, $4 - 1 + $9, strand, umi_match[1];
  }
}
