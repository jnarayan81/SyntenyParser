# SyntenyParser

maf2syntenyparser script usage

```
  print "Usage: $0 --qfile --flip 1 \n\n";
  print	"Options:\n";
  print "	--qfile|-q	query Synteny blocks_coords file\n";
  print "	--flip|-f	flip the coordiantes if negative oriented | 1 for yes or 0 for no \n";
  print "	--core|-c	Number of core/core to use \n";
  print "	--ofile|-o	oufile/results file\n";
  print "	--ref|-r	reference name \n";
  print "	--tar|-t	target name\n";
  print "	--plot|-p	name of scaff to plot\n";
  print "	--block|-b	report block format : for default see below\n";
  print "     	--help|-h	brief help message\n";

```


Outfile format (default)
```
  +Col1		Name of the reference
  +Col2		Name of the reference chromosome/contig/sequence
  +Col3		Strand
  +Col4		Start coordinates
  +Col5		End coordinates
  +Col6		Length of alignments
  +Col7		Name of the target
  +Col8		Name of the target chromosome/contig/sequence
  +Col9		Strand
  +Col10	Start coordinates
  +Col11	End coordinates
  +Col12	Length of alignments
  +Col13	Others

```

