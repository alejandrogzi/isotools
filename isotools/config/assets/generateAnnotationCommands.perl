#!/usr/bin/perl
# Michael Hiller, 2024
# 

use strict;
use warnings;
use diagnostics;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Scalar::Util::Numeric qw(isint);

$| = 1;		# == fflush(stdout)
my $verbose = 0;

# Globals:
my %defVars = ();       # hash with variables defined in DEF

my $genomeFasta = "";   # pointing to the genome fasta file: $workDir/genome.fa
my $genome2bit = "";    # pointing to the genome 2bit file: $workDir/genome.2bit
my $assembly = "";      # empty, unless 'genome' parameter was the name to an assembly (e.g. hg38)

# which steps to omit
my $skipCCS = 0;
my $skipLima = 0;
my $forceOverwrite = 0;

my $call = "";   # a command to be executed via system() or backticks
my $genomePath = $ENV{'genomePath'};   # should be set on the hillerlab system. Only required when genome is an assembly name

# output script containing all commands
my $outputScript = "doAnnotation.sh";

# options
my $usage = "usage: $0 genome DEF outputDir [Optional parameters]\n
Generates a shell script '$outputScript' that runs all steps.
$0 itself does not execute any of these steps. 
$0 must be executed on delta and in /genome/gbdb-HL/$assembly/GenomeBuild !! 
\
genome                      3 options are possible
                            1. /path/to/genome.fa file. File must be unzipped
			    2. /path/to/genome.2bit file. twoBitToFa will be used to convert to a .fa file
			    3. assembly name. In this case there should be a \$genomePath/gbdb-HL/\$genome/\$genome.2bit file. 
DEF                         DEF file listing all parameters and options 
outputDir                   path or name of output directory (e.g. annotationResults)
                            script will error-exit if that directory already exists. Use -f to overwrite previous results
\
Optional parameters
-skipCCS                    flag: If set, do not run CCS
-skipLima                   flag: If set, do not run Lima
-f                          force overwrite file in outputDir if that exists
-v                          verbose output
\n";

GetOptions ("v|verbose"  => \$verbose, "skipCCS" => \$skipCCS, "skipLima" => \$skipLima, "f" => \$forceOverwrite) || die "$usage\n";
die "$usage\n" if ($#ARGV < 2);

my $genome = $ARGV[0];  # can be a fasta file (.fa), a 2bit file (.2bit) or an assembly name
my $DEF = $ARGV[1];
my $outDir = $ARGV[2];


loadDef($DEF);
checkDef();

# check whether we are on delta
die "ERROR: $0 must be executed on delta\n" if ($ENV{'HOSTNAME'} ne "delta");
my $curDir = `pwd`; chomp($curDir);

# outDir
die "ERROR: $outDir already exists. Either set aside or delete. Or use -f to overwrite\n" if ((-d $outDir) && ($forceOverwrite == 0));
$call = "mkdir -p $outDir";
system("$call") == 0 || die "ERROR: $call failed\n";

# working directory
$call = "date \"+%Y_%m_%d_%H.%M\"";
my $date=`$call`; die "ERROR: $call FAILED\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);  
chomp($date);
$call = "mktemp -d $defVars{'TMPDIR'}/TEMPAnnotation_${date}_XXX";
my $workDir = `$call`; die "ERROR: $call FAILED\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
chomp($workDir);
print "created temp dir: $workDir\n";


# check or create genome fasta file
# option 1: user parameter is a 2bit file
if ($genome =~ /.2bit$/) {
	die "ERROR: $genome appears to be a 2bit file (ending .2bit) but the file does not exist\n" if (! -f $genome);
	convertTwoBitToFa($genome, "$workDir/genome.fa");
	$genomeFasta = "$workDir/genome.fa";
	# link 2bit
	my $fullGenomeDir = getFullDir($genome);   # return full path if a relative path is given
	$call = "ln -s $fullGenomeDir $workDir/genome.2bit";
	system("$call") == 0 || die "ERROR: $call failed\n";
	$genome2bit = "$workDir/genome.2bit";
# option 2: user parameter is a fasta file
}elsif ($genome =~ /.+\.fa*/) {
	$genomeFasta = $genome;
	# must exist and be unzipped
	die "ERROR: $genomeFasta does not exist\n" if (! -f $genomeFasta);
	die "ERROR: $genomeFasta must be an unzipped file\n" if ($genomeFasta =~ /\w+\.gz$/);
	die "ERROR: $genomeFasta is not a fasta file\n" if (`cat $genomeFasta | head -n1` !~ /^>/); 
	convertFaToTwoBit($genome, "$workDir/genome.2bit");
	$genome2bit = "$workDir/genome.2bit";
	# link 2bit
	my $fullGenomeDir = getFullDir($genome);   # return full path if a relative path is given
	$call = "ln -s $fullGenomeDir $workDir/genome.fa";
	system("$call") == 0 || die "ERROR: $call failed\n";
	$genomeFasta = "$workDir/genome.fa";
# option 3: user parameter is an assembly name (e.g. hg38)
}else {
	die "ERROR: genome ($genome) appears to be an assembly name (no .2bit or .fa suffix), but environment variable genomePath is not set\n" if ($genomePath eq "");
	$genome2bit = "$genomePath/gbdb-HL/$genome/$genome.2bit";
	die "ERROR: cannot access assembly 2bit at $genome2bit\n" if (! -f $genome2bit);
	convertTwoBitToFa($genome2bit, "$workDir/genome.fa");
	$genomeFasta = "$workDir/genome.fa";
	# link 2bit
	$call = "ln -s $genome2bit $workDir/genome.2bit";
	system("$call") == 0 || die "ERROR: $call failed\n";
	$genome2bit = "$workDir/genome.2bit";
	$assembly = $genome;    # assembly now points to an assembly name. Can be used for loading tracks into the browser
}


my $baseDir = dirname($assembly);
my $assemblyName = basename($assembly);

#die "ERROR: $ENV{'HOME'}/src/userBrowserTracks/hillerlab/$trackDbDirectory/ does not exist\n" if (! -d "$ENV{'HOME'}/src/userBrowserTracks/hillerlab/$trackDbDirectory/");
print "output directory: $outDir   assembly: $assemblyName\n";


# get email of the user to send an email at the end
$call="getent passwd \"\$USER\" | cut -f5 -d \":\" | sed 's/ /./g'";
my $username=`$call`;
die "ERROR: $call FAILED\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
chomp($username);
my $email="$username\@senckenberg.de";
print "your email is: $email\n\n";



# start generating the script
open(outFile, ">$outputScript") || die "ERROR: cannot write to $outputScript\n";
print outFile "#!/bin/bash\nset -eux\nset -o pipefail\n\n";

print outFile "trap \"mailx -a /genome/gbdb-HL/$assembly/GenomeBuild/log.txt -s \\\"$outputScript $assembly CRASHED. Log is attached\\\" $email < /dev/null\" ERR\n";



# raw data; Use CCS to convert the subreads into HiFi reads (has adapters and polyA)
if (! $skipCCS) {
	print outFile "module load pbccs/$defVars{'pbccs'}\n";
	print outFile "ccs\n";
}

# Use Lima to demultiplex barcodes and remove primers to produce full-length (FL) inserts (NOTE: full length inserts, not necessarily FL transcripts)
if (! $skipLima) {
	print outFile "module load lima/$defVars{'lima'}\n";
	print outFile "lima\n";
}


# isoseq3 refine
print outFile "module load isoseq/$defVars{'isoseq'}\n";

=pod
ccs (pbccs v.4.2.0-1)	ccs ${raWData_Link}/${dataPrefix}.subreads.bam ${P_out_CCS}/${dataPrefix}.ccs.$i.bam --min-rq 0.95 --chunk $i/$nChunks -j $defVars{'nThreadsIsoSeq3'} --report-file ${runID}_ccs_report.txt
lima v.1.11.0	lima --isoseq --peek-guess ${P_out_CCS}/${dataPrefix}.ccs.${i}.bam ${f_primers} ${P_out_lima}/${dataPrefix}.fl.${i}.bam
isoseq3 refine v.3.4.0	isoseq3 refine ${P_out_lima}/${dataPrefix}.fl.${i}.${primer_p5}${primer_p3}.bam ${f_primers} ${P_out_isoR}/${dataPrefix}.flnc.${i}.bam --num-threads $defVars{'nThreadsIsoSeq3'}

# cluster. Important, we keep low coverage transcripts (singletons), as this pipeline has an advanced transcript classification procedure.
# : Hierarchical n*log(n) clustering, meaning iterative alignment of shorter to longer 
isoseq3 cluster v.3.4.0	isoseq3 cluster ${P_out_isoC}/ALL.flnc.fofn ${P_out_isoC}/ALL.CuP.bam --singletons -verbose --split-bam $splitBAM --num-threads $defVars{'nThreadsIsoSeq3'} --log-file ${P_out_isoC}/ALL.CuP.log

# minimap2 to map the reads against the assembly
my $maxIntronSize = "1600k";

#print outFile "minimap2 --eqx -a -c -t $defVars{'nThreadsMinimap'} -ax splice:hq -uf --secondary=no -C5 -o ${P_outMinimap2}/ALL.CuP.aln.sam --junc-bed ${ref_annot} -cs long -G $maxIntronSize $genomeFasta ${P_out_isoC}/ALL.CuP.fasta.gz\n";

# first run the filter script but have more lax parameters for %id and 3'clip
print outFile "filterMinimapQuality.perl $input.sam -perID 96 -clip3 50 -polyAReadSuffix 30\n";
print outFile "correctMinimap $pathToTOGAQueryAnnotation.bed $input.good.sam $assembly $input.correctedMinimap.sam\n";
# now with default %id 98 and max 3'clip 20
print outFile "filterMinimapQuality.perl $input.correctedMinimap.sam -polyAReadSuffix 30\n";
# resulting file to move on with is called $input.correctedMinimap.good.sam
print outFile "filterArtificialLongIntrons ....\n";

=cut


=pod
# load
my $gapParameter = "";
$gapParameter = "-g $gapFile" if ($gapFile ne "");
# Note that we tested above that gapFile and minAssemblyGapSize are not both specified
$gapParameter = "-M $minAssemblyGapSize" if ($minAssemblyGapSize != -1 && $minAssemblyGapSize != $defaultminAssemblyGapSize); 
print outFile "loadGenome.sh -i config.ra -f $genomeFasta -m $assembly.sorted.fa.out.gz $gapParameter\n";
print outFile "djStuff/cleanUp_$assembly.sh\n\n";

# load the optional telomere file if parameter is given
if ($telomereFile ne "") {
        print outFile "hgLoadBed $assembly HLtelomere $telomereFile\n";
}

# alignment
if ($doAlignment ne "") {
	my @refs = split(/,/, $doAlignment);
	foreach my $ref (@refs) {
		print "output alignment commands for reference $ref\n" if ($verbose);
		print outFile "ssh delta \"makePairwiseAlignments.sh -r $ref -a $assembly\"\n"; 
		print outFile "\n";
	}

	foreach my $ref (@refs) {
		print outFile "cd /genome/gbdb-HL/$ref/lastz/\n";
		print outFile "rsync -av delta:/projects/hillerlab/genome/gbdb-HL/$ref/lastz/vs_$assembly/ vs_$assembly/\n";
		print outFile "cd vs_$assembly/axtChain\n";
		print outFile "hgLoadChain -tIndex $ref HLchainFilled$assembly $ref.$assembly.allfilled.chain.gz\n";
		print outFile "hgLoadNet $ref HLnetFilled$assembly $ref.$assembly.net.gz\n";
		print outFile "\n";
	}
}

# alignment
if ($doTOGA ne "") {
	my @refs = split(/,/, $doTOGA);
	foreach my $ref (@refs) {
		print "output TOGA commands for reference $ref\n" if ($verbose);
		print outFile "ssh delta \"do_all_togaHL.py $ref $assembly\"\n";
	}
	print outFile "\n";
	print outFile "echo \"check if the TOGA bigBed files are placed in the /var/www/data directory\"\n";
	print outFile "ls -l /var/www/data/$assembly\n";
	print outFile "\n";
}

# create userBrowserTracks folder and trackDb.ra
print outFile "cd ~/src/userBrowserTracks/hillerlab/\n";
print outFile "git pull\n";
print outFile "mkdir $trackDbDirectory/$assembly\n";
print outFile "touch $trackDbDirectory/$assembly/trackDb.ra\n";
print outFile "git add $trackDbDirectory/$assembly/trackDb.ra\n";
print outFile "git commit -m \"added $assembly\"\n"; 
print outFile "git push\n";
print outFile "make $assembly\n";
print outFile "\n";

# update hg38 / mm10 / galGal6 chainNet master track
if ($doAlignment ne "") {
	my @refs = split(/,/, $doAlignment);
	foreach my $ref (@refs) {
		if ($doAlignment =~ /hg38/ || $doAlignment =~ /mm10/) {
			# call create_ChainNet_Mammals.py for hg38 or mm10 
			# --> for other mammals, we don't have a composite track setup right now
			print outFile "cd ~/src/userBrowserTracks/hillerlab/\n";
			print outFile "src/create_ChainNet_Mammals.py -r $ref\n";
			print outFile "./$ref.chainNet.mammals.README\n";
			print outFile "\n";
		}elsif ($doAlignment =~ /galGal6/) {
			# call create_ChainNet_Birds.py for galGal6 
			# --> for other birds, we don't have a composite track setup right now
			print outFile "cd ~/src/userBrowserTracks/hillerlab/\n";
			print outFile "src/create_ChainNet_Birds.py -r $ref\n";
			print outFile "./$ref.chainNet.README\n";
			print outFile "\n";		
		}
	}


}

=cut
print outFile "\necho \"ALL SUCCESSFULLY DONE\"\n";

print outFile "mailx -a /genome/gbdb-HL/$assembly/GenomeBuild/log.txt -s \"$outputScript $assembly successfully finished. Log is attached\" $email < /dev/null\n";

close outFile;

`chmod +x $outputScript`;
print "\n==> Successfully generated $outputScript\nPlease carefully check this file. If everything is correct, run\nscreen -S $assembly; $outputScript >> log.txt 2>> log.txt\n";





# wrapper around twoBitToFa 
sub convertTwoBitToFa {
	my ($twoBitFile, $outFa) = @_;
	my $call = "twoBitToFa $twoBitFile $outFa"; 
	system("$call") == 0 || die "ERROR: $call failed\n";
}

# wrapper around faToTwoBit 
sub convertFaToTwoBit {
	my ($faFile, $twoBitFile) = @_;
	my $call = "faToTwoBit $faFile $twoBitFile"; 
	system("$call") == 0 || die "ERROR: $call failed\n";
}


#########################################################################
# The following routines were taken almost verbatim from UCSCs blastz-run-ucsc
#########################################################################
sub loadDef {
  # Read parameters from a bash script
  my ($def) = @_;
  my $fh = mustOpen("$def");
  while (<$fh>) {
    s/^\s*export\s+//;
    next if (/^\s*#/ || /^\s*$/);
    if (/(\w+)\s*=\s*(.*)/) {
      my ($var, $val) = ($1, $2);
      while ($val =~ /\$(\w+)/) {
	my $subst = $defVars{$1};
	if (defined $subst) {
	  $val =~ s/\$$1/$subst/;
	} else {
	  die "Can't find value to substitute for \$$1 in $DEF var $var.\n";
	}
      }
      $defVars{$var} = $val;
    }
  }
  close($fh);
  
  # test if TMPDIR in DEF exists; if not create it
#  if (exists $defVars{TMPDIR} && ! -e $defVars{TMPDIR}) {
#  		carp "create $defVars{TMPDIR}\n";
#		system "set -o pipefail; mkdir -p $defVars{TMPDIR}";
#  }
  
}

sub mustOpen {
  # Open a file or else die with informative error message.
  my ($fileSpec) = @_;
  die "Must have exactly 1 argument" if (scalar(@_) != 1);
  die "undef input" if (! defined $fileSpec);
  open(my $handle, $fileSpec)
    || die "Couldn't open \"$fileSpec\": $!\n";
  return $handle;
}


sub checkDef {
  # Make sure %defVars contains what we need and looks consistent with
  # our assumptions.
  foreach my $req ('TMPDIR', 'nThreadsIsoSeq3', 'nThreadsMinimap') {
      &requireVar("$req");
  }
#    &requirePath($s . 'DIR');

#  $tDb = &getDbFromPath('SEQ1_DIR');
#  $qDb = &getDbFromPath('SEQ2_DIR');

  print "$DEF looks OK!\n";
}

sub requireVar {
  my ($var) = @_;
  die "Error: $DEF is missing variable $var\n" if (! defined $defVars{$var});
}

sub getFullDir {
  my $dir = shift;
  $call = "realpath $dir";
  return $dir if ($dir =~ /^\//);  # is already an absolute path
  # otherwise call Linux realpath
  my $fullDir = `$call`;
  die "ERROR: $call FAILED. Maybe the 'realpath' linux command does not exist on your system? If so, provide a full path for $dir\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0); 
  chomp($fullDir);
  return $fullDir;
}  

