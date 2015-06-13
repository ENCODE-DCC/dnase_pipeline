table narrowPeak
"ENCODE narrowPeak"
(
string  chrom;          "Reference sequence chromosome or scaffold"
uint    chromStart;     "Start position of feature on chromosome"
uint    chromEnd;       "End position of feature on chromosome"
string  name;           "Name of gene"
uint    score;          "Score"
char[1] strand;         "+ or - for strand"
float   SignalValue;     "Measurement of overall (usually, average) enrichment for the region. "
float  pValue;       "-log10 p value"
float  qValue;       "-log10 q value"
int    peak;         "Point-source called for this peak"
)