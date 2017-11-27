#!/usr/bin/perl

use strict;
use 5.010;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Pod::Usage;
use Parallel::ForkManager;
use Chromosome::Map;
 
# maf2synteny parser
#Author: Jitendra Narayan
#Usage: perl SyntenyParser.pl -a blocks_coords.txt -f 1 -c 40 -o sampleOut -r ref -t tar -p scaffold_2

#print "\nNOTE: The script assume each blocks number has two only lines in your blocks_coords.txt file\n";

my ($afile, $help, $flip, $man, $core, $ofile, $ref, $tar, $block, $plot, $mode, $chrlen);
my $version=0.1;
GetOptions(
    'afile|a=s' => \$afile,
    'flip|f=n' => \$flip,
    'ref|r=s' => \$ref,
    'tar|t=s' => \$tar,
    'block|b=s' => \$block,
    'plot|p=s' => \$plot,
    'core|c=n' => \$core,
    'ofile|o=s' => \$ofile,
    'mode|m=s' => \$mode,
    'lfile|l=s' => \$chrlen,
    'help|h' => \$help
) or die &help($version);

&help($version) if $help;
#pod2usage("$0: \nI am afraid, no files given.")  if ((@ARGV == 0) && (-t STDIN));
my $coreNum = `grep -c -P '^processor\\s+:' /proc/cpuinfo`;
if ($core > $coreNum)  {print "Arrr ! Are you kidding me, you have only core $coreNum\n"; exit;}
if (!$afile or !$flip or !$core or !$ofile or !$mode) {help($version) }

if ($core > 1) { print "\nWARNING: Might distort the outfile, use more core/thread at your own risk\n---Advised to use 1 thread for now---\n";}

if ($block and $plot) { print "\nSorry plotting function is available with default option of blocks creation\n Try without -b flag\n"; exit; }


if (lc($mode) eq "sibelia") { 
#Seperator
local $/ = "--------------------------------------------------------------------------------";

open my $afh, '<', $afile or die "Unable to open $afile: $!\n";
open my $ofh, '>', $ofile or die "Unable to open $ofile: $!\n";
my %blocks; my %terms; my $blkCnt;
while (<$afh>) {
    chomp;
    next if /^\s*$/;
    #Lets store the ids
    if ( $. == 1 ) { %terms = storeIds($_); next;}
    $blkCnt++;
    $blocks{$blkCnt} = $_;
    #printBlocks ($_, \%terms);
}
close $afh;

print "\nMaf2Synteny Parsing file finished\nTotal number of blocks $blkCnt\nWorking for formating ...\n";

my $max_procs = $core;
my @blkCntNo = keys %blocks;
# hash to resolve PID's back to child specific information
my $pm =  new Parallel::ForkManager($max_procs);
# Setup a callback for when a child finishes up so we can
# get it's exit code
  $pm->run_on_finish (
    sub { my ($pid, $exit_code, $ident) = @_;
      #print "** $ident just got out of the pool ". "with PID $pid and exit code: $exit_code\n";
    }
  );

  $pm->run_on_start(
    sub { my ($pid,$ident)=@_;
     #print "** $ident started, pid: $pid\n";
    }
  );

  $pm->run_on_wait(
    sub {
      #print "** Have to wait for one children ...\n"
    },
    0.5
  );

  NAMES:
  foreach my $child ( 0 .. $#blkCntNo ) {
    #next if length($blocks{$blkCntNo[$child]}) <= $length;
    my $pid = $pm->start($blkCntNo[$child]) and next NAMES;
    printBlocks($blocks{$blkCntNo[$child]}, \%terms);
    $pm->finish($child); # pass an exit code to finish
  }
  print "Waiting for all the jobs to complete...\n";
  $pm->wait_all_children;
  print "DONE ... Everybody is out of the computation pool!\n";
  close $ofh;


sub printBlocks {
my ($line, $terms_ref)= @_;
my %terms=%$terms_ref;
my @val = split /\n/, $line;
shift @val for 1..3; # Delete the "Blocks# and "Header"
my $refBlk = shift @val;
my @refBlkVal = split /\t/, $refBlk;
if (!$ref) {$ref="REF";}
if (!$tar) {$tar="TAR";}
my $refLine='NA';

if ($block eq 'satsuma') { $refLine="$terms{$refBlkVal[0]}\t$refBlkVal[2]\t$refBlkVal[3]"; } 
elsif ($block eq 'hsb') { $refLine="$ref\t$terms{$refBlkVal[0]}\t$refBlkVal[2]\t$refBlkVal[3]"; }
else { $refLine="$ref\t$terms{$refBlkVal[0]}\t$refBlkVal[1]\t$refBlkVal[2]\t$refBlkVal[3]\t$refBlkVal[4]"; }
my $refC = $terms{$refBlkVal[0]};
my $sname='targetSpsName';
foreach (@val) {
	my @v = split /\t/;
	my @nan= split /\:/, $terms{$v[0]};
	my $st="NA"; my $ed="NA";
	if (($flip) and ($v[1] eq "-")) { $st=$v[3]; $ed=$v[2];} else { $st=$v[2]; $ed=$v[3];}
	my $selfDec='aln';
	if ($refC eq $terms{$v[0]}) { $selfDec='self';}
	#satsuma#scaffold_1_1087316_bp	1	2945	scaffold_1-edited	1	2937	0.984715	+
	if ($block eq 'satsuma') { print $ofh "$refLine\t$terms{$v[0]}\t$st\t$ed\t$v[4]\t$v[1]\n"; }
	#hsb#galgal:100k	1	31014	8882608	1	5187	8577687	+	anas_platyrhynchos	chromosomes
	elsif ($block eq 'hsb') {  print $ofh "$refLine\t$tar\t$terms{$v[0]}\t$st\t$ed\t$v[1]\t$sname\tchr/scaff\n"; }
	#default format see README
	else { print $ofh "$refLine\t$tar\t$terms{$v[0]}\t$v[1]\t$st\t$ed\t$v[4]\t$selfDec\n"; }
	}
}
}

elsif (lc($mode) eq "lastz") {

#check if user have provided length file
if (!$chrlen) { print "It seems you forgot to provide length file\n"; }
#remove overlapping aln
uniqAln($afile, 'tmpAln');
my $lenHash_ref=storeLen($chrlen);
my %lenHash = %$lenHash_ref;
#print the final format
open my $afh, '<', 'tmpAln' or die "Unable to open tmpAln: $!\n";
open my $ofh, '>', $ofile or die "Unable to open $ofile: $!\n";
while (<$afh>) {
    chomp;
    next if /^\s*$/;
    #Lets store the ids
    my @aLine = split '\t', $_;
    my $refL=$aLine[5]-$aLine[4];
    my $tarL=$aLine[10]-$aLine[9];
    my $refLen=$lenHash{$aLine[1]};
    print $ofh "$ref\t$refLen:$aLine[1]\t$aLine[2]\t$aLine[4]\t$aLine[5]\t$refL\t$tar\t$aLine[6]\t$aLine[7]\t$aLine[9]\t$aLine[10]\t$tarL\taln\n";
}
close $afh;
}


#Lets begin drawing
print "Lets plot the graph for chr/contig/scaff $plot : \n";
 my ($sHash_ref, $cSize) = extract2Map($ofile,$plot);
#my $cSizeKB=$cSize/100;
 plotSyn('tmp.aln', $plot, $cSize, $sHash_ref);





#subs here --------------------------------------------------------------------
#store len file
sub storeLen {
my ($lFile)=@_;
my %lHash;
open my $lfh, '<', $lFile;
while (<$lfh>) {
    chomp;
    next if /^\s*$/;
    my @tmpLen = split '\t', $_;
    $lHash{$tmpLen[0]}=$tmpLen[1];
}
return \%lHash;
}

#create uniq aln file
sub uniqAln {
my ($iName, $oName)=@_;
open my $ifh, '<', $iName;
open my $ofh, '>', $oName;

my @terms;
while (<$ifh>) {
    chomp;
    push @terms, [split /\t/];
}

my $biggest = 0;
my $id = '';

for my $term (sort sorter @terms) {
	$biggest = 0 if $id ne $term->[1];
    if ($term->[5] > $biggest) {
        my $nLine= join "\t", @$term;
	print $ofh "$nLine\n";
        $biggest = $term->[5];
    }
    $id = $term->[1];    
}
}

#sort the aln
sub sorter {
	$a->[1] cmp $b->[1] ||
     $a->[4] <=> $b->[4]
  || $b->[5] <=> $a->[5]
}

#plot synteny
sub plotSyn {
my ($tmpFile, $name, $len, $sHash_ref)=@_;
 my %sHash=%$sHash_ref;
my $map = Chromosome::Map->new (-length     => $len, 
                                 -name       => $name, 
                                 -height     => '1500', 
                                 -units      => 'bp', 
                                ); 
 
my $size  = $map->get_map_size; 
my $units = $map->get_map_units; 
print "Map size: $size $units\n";
 
my $qtl_track  = Chromosome::Map::Track->new (-name => 'QTL',
                                               -type => 'interval',
                                              );
# adding tracks to map
$map->add_track($qtl_track);
 
my $nb_track = $map->get_nb_tracks;
print "Nb track: $nb_track\n";
 
# my @Color = qw (blueviolet darkgoldenrod black softblue khaki red blue tomato);

#plot all blocks against chromosome/contigs/scaffolds
foreach my $sBlocks (keys %sHash) {
     my @sBlocks = split '\t', $sHash{$sBlocks};
     my @tmpName=split '\:', $sBlocks[7]; #Tar name
     my $localName="$tmpName[1]:$sBlocks[3]-$sBlocks[4]";
     my $lColor='black';
     if ($sBlocks[8] eq '-') { $lColor= 'red';} else {$lColor= 'darkgoldenrod';}
     my $qtl1 = Chromosome::Map::Block->new (-name  => $localName,
                                         -start => $sBlocks[3],
                                         -end   => $sBlocks[4],
                                         -color => $lColor,
                                        );
     $qtl_track->add_element($qtl1);
}
 
my $png = $map->png;
my $filename_png = "chr_map_$name.png";
open (PNG, ">$filename_png") || die "cannot create file: $filename_png!\n";
binmode PNG;
print PNG $png;
close PNG;

}

#Lets plot the map
sub extract2Map {
my ($alnFile, $plot) = @_;
my %subHash; my $cLen;
my $tmpfile='tmp.aln';
local $/ = "\n";
open my $alnfh, '<', $alnFile or die "Unable to open $alnFile: $!\n";
open my $alnofh, '>', $tmpfile or die "Unable to open $tmpfile: $!\n";
while (<$alnfh>) {
    chomp;
    next if /^\s*$/;
    my @aln= split /\t/;
    my @rName = split '\:', $aln[1];
    next if $aln[12] eq 'self';
    if ($plot eq $rName[1]) {
	print $alnofh "$_\n";
	$subHash{$.}=$_;
        $cLen=$rName[0];
	}
}
if (!$cLen) { print "Length for $plot is empty -- Check you requested chr/scaff/contig name\n"; exit;}
return (\%subHash, $cLen);
}

#store subs
sub storeIds {
my $line = shift;
my %terms;
my @val = split /\n/, $line;
foreach (@val) { my @v = split /\t/; $terms{$v[0]} = "$v[1]:$v[2]";} 
return %terms;
}

#Help section
sub help {
  my $ver = $_[0];
  print "\n SyntenyParser.pl $ver\n";
  print "\n Report errors/bug to Jitendra 'jnarayan81ATgmail.com'\n\n";

  print "Usage: $0 --afile --flip 1 \n\n";
  print	"Options:\n";
  print "	--afile|-a	alignment 'Maf2Synteny blocks_coords/lastz general' file\n";
  print "	--flip|-f	flip the coordiantes if negative oriented | 1 for yes or 0 for no \n";
  print "	--core|-c	Number of core/core to use \n";
  print "	--ofile|-o	oufile/results file\n";
  print "	--ref|-r	reference name \n";
  print "	--tar|-t	target name\n";
  print "	--plot|-p	name of scaff to plot\n";
  print "	--block|-b	report block format\n";
  print "	--lfile|-l	provide chr length file\n";
  print "	--mode|-m	provide alingnment mode\n";
  print "     	--help|-h	brief help message\n";

print "For LastZ general-: perl SyntenyParser.pl -a seeALN_scaffold_15.lz -f 1 -c 1 -o see2 -r ref -t tar -p scaffold_15 -m lastz -l scaff15.fa.fai
";
print "For Sibelia: perl SyntenyParser.pl -a blocks_coords.txt -f 1 -c 1 -o see2 -r ref -t tar -p scaffold_1 -m sibelia";
exit;
}

