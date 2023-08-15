#!/usr/bin/perl

BEGIN {
	use FindBin '$RealBin';
	push @INC, $RealBin;
}

use Getopt::Std; %opt = ();

getopts("O:t:sm:G:C:", \%opt);

# dfaults
$minSize = 10000000;

$die = "

sciMET_BSBolt2cellCalls.pl (options) [rmdup & filtered bam file, name-sorted]

Extracts methylation from a BSBolt bam and outputs bismark annotated
calls for each cell.

Options:
   -O   [STR]   Out prefix
   -m   [INT]   Minimum chromosome size to retain (def = $minSize)
                  Used to exclude random and other small contigs
   -t   [INT]   Max number of concurrent threads (def = 1)
   -C   [STR]   List of cellIDs to include (use instead of pre-
                  filtering bam file)
   
";

if (!defined $ARGV[0]) {die $die};

if (defined $opt{'C'}) {
	open IN, "$opt{'C'}";
	while ($cellID = <IN>) {
		chomp $cellID;
		$INCLUDE{$cellID} = 1;
	} close IN;
}

%CALL_CONV = ("x" => "z", "X" => "Z", "y" => "x", "Y" => "X", "z" => "h", "Z" => "H");

if (!defined $opt{'s'}) { # main thread
	#defaults
	if (!defined $opt{'t'}) {$opt{'t'} = 1};
	
	# open log file
	open LOG, ">$opt{'O'}.log";
	$ts = localtime(time);
	print LOG "$ts\tProgram called.\n\n============== Phase 1: Parsing bam and processing per-cell calls ==============\n\n";

	# setup output directory
	system("mkdir $opt{'O'}");
	$ts = localtime(time);
	print LOG "$ts\tOutput directory created: $opt{'O'}\n";
	
	# setup threading stuff
	%QUEUE = ();
	
	# start parsing bam
	$ts = localtime(time);
	print LOG "$ts\tParsing bam file: $ARGV[0]\n";
	
	open IN, "samtools view $ARGV[0] |";
	$currentBarc = "null";
	$methCol = "null";
	
	$prevID = "null";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$barc = $P[0]; $barc =~ s/:.+$//;
		$readID = $P[0]; $readID =~ s/#.+$//;
		
		# check if barc is included, ortherwise skip line
		if (!defined $opt{'C'} || defined $INCLUDE{$barc}) {
		
			# check to close out barcode
			if ($currentBarc ne $barc) {
				# close out prev barc
				if ($currentBarc ne "null") {
					close OUT;
					$QUEUE{$currentBarc} = 0;
				}
				# check thread stats
				check_threads();
				
				# set next barcode & open outfile
				$currentBarc = $barc;
				$prevID = "null";
				open OUT, ">$opt{'O'}/$currentBarc.meth";
			}
			# parse read
			$chr = $P[2];
			$pos = $P[3];
			# ID column in bam that has the meth col if not yet found
			if ($methCol eq "null") {
				for ($i = 11; $i < @P; $i++) {
					if ($P[$i] =~ /^XB/) {
						$methCol = $i;
						$i+=100;
					}
				}
			}
			$meth = $P[$methCol];
			$meth =~ s/^XB:Z://;
			
			if ($prevID eq "null") {
				$prevChr = $chr;
				$prevPos = $pos;
				$prevMeth = $meth;
				$prevID = $readID;
			} elsif ($prevID eq $readID) { # print pair
				print OUT "$prevChr\t$prevPos\t$prevMeth\t$chr\t$pos\t$meth\n";
				$prevID = "null";
			} else { # solo - print previous as solo & save new
				print OUT "$prevChr\t$prevPos\t$prevMeth\n";
				$prevChr = $chr;
				$prevPos = $pos;
				$prevMeth = $meth;
				$prevID = $readID;
			}
			
			#print OUT "$chr\t$pos\t$meth\n";
			
		}
	} close IN;
	
	# finish last one
	close OUT;
	$QUEUE{$currentBarc} = 0;
	
	# process rest of queue
	$incomplete = 1;
	while ($incomplete > 0) {
		$incomplete = 0;
		sleep(10); # wait
		check_threads();
		foreach $queued_barc (keys %QUEUE) {
			if ($QUEUE{$queued_barc} < 2) { # active or waiting
				$incomplete++;
			}
		}
	}
	$ts = localtime(time);
	print LOG "$ts\tAll threads complete! DONE!\n";


} else { # subthread
	open IN, "$ARGV[0]/$ARGV[1].meth";
	open CG, "| sort -k1 -k 2,2n > $ARGV[0]/$ARGV[1].CG.cov";
	open CH, "| sort -k1 -k 2,2n > $ARGV[0]/$ARGV[1].CH.cov";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\s/, $l);
		
		%read_cov = (); %read_meth = ();
		if (defined $P[3]) { # pair
			load_read($P[0],$P[1],$P[2]);
			load_read($P[3],$P[4],$P[5]);
		} else {
			load_read($P[0],$P[1],$P[2]);
		}
		
		foreach $coord (keys %read_cov) {
			if (!defined $COORD_cov{$coord}) {
				if (defined $read_cov{$coord}{'x'}) {$COORD_cov{$coord}{'x'} = 1};
				if (defined $read_cov{$coord}{'h'}) {$COORD_cov{$coord}{'h'} = 1};
				if (defined $read_meth{$coord}{'x'}) {$COORD_meth{$coord}{'x'} = 1};
				if (defined $read_meth{$coord}{'h'}) {$COORD_meth{$coord}{'h'} = 1};
			} else {
				if (defined $read_cov{$coord}{'x'}) {$COORD_cov{$coord}{'x'}++};
				if (defined $read_cov{$coord}{'h'}) {$COORD_cov{$coord}{'h'}++};
				if (defined $read_meth{$coord}{'x'}) {$COORD_meth{$coord}{'x'}++};
				if (defined $read_meth{$coord}{'h'}) {$COORD_meth{$coord}{'h'}++};
			}
		}
		
	} close IN;
	
	foreach $coord (keys %COORD_cov) {
		if (defined $COORD_cov{$coord}{'x'}) {
			if (!defined $COORD_meth{$coord}{'x'}) {$COORD_meth{$coord}{'x'}=0};
			$pct = sprintf("%.2f", ($COORD_meth{$coord}{'x'}/$COORD_cov{$coord}{'x'})*100);
			$unmeth = $COORD_cov{$coord}{'x'}-$COORD_meth{$coord}{'x'};
			print CG "$coord\t$pct\t$unmeth\t$COORD_meth{$coord}{'x'}\n";
		}
		if (defined $COORD_cov{$coord}{'h'}) {
			if (!defined $COORD_meth{$coord}{'h'}) {$COORD_meth{$coord}{'h'}=0};
			$pct = sprintf("%.2f", ($COORD_meth{$coord}{'h'}/$COORD_cov{$coord}{'h'})*100);
			$unmeth = $COORD_cov{$coord}{'h'}-$COORD_meth{$coord}{'h'};
			print CH "$coord\t$pct\t$unmeth\t$COORD_meth{$coord}{'h'}\n";
		}
	}
	
	$ts = localtime(time);
	system("echo '$ts' > $ARGV[0]/$ARGV[1].complete");
	exit;
}

sub load_read {
	$chr = $_[0]; $pos = $_[1]; $meth = $_[2];
	@M = split(//, $meth);
	for ($i = 0; $i < @M; $i++) {
		if ($M[$i] !~ /[0-9]/) {
			
			$coord = "$chr\t$pos";
			if ($M[$i] eq "x") { #CG, unmeth
				$read_cov{$coord}{'x'} = 1;
			} elsif ($M[$i] eq "X") { #CG, meth
				$read_cov{$coord}{'x'} = 1;
				$read_meth{$coord}{'x'} = 1;
			} elsif ($M[$i] eq "y" || $M[$i] eq "z") { # CH , unmeth
				$read_cov{$coord}{'h'} = 1;
			} elsif ($M[$i] eq "Y" || $M[$i] eq "Z") { # CH , meth
				$read_cov{$coord}{'h'} = 1;
				$read_meth{$coord}{'h'} = 1;
			}
			
			$pos++;
			
		} elsif ($M[$i] =~ /[0-9]/) {
			$add = $M[$i];
			while ($M[$i+1] !~ /xyz/i && $M[$i+1] =~ /[0-9]/) {
				$i++;
				$add .= $M[$i];
			}
			$pos += $add;
		}
	}
}
	
sub check_threads {
	$running = 0; @WAITING = ();
	foreach $queued_barc (keys %QUEUE) {
		if ($QUEUE{$queued_barc} == 1) { # active
			$running++;
			if (-e "$opt{'O'}/$queued_barc.complete") { # finished
				$running--;
				$QUEUE{$queued_barc} = 2;
				print LOG "\t\t\t$queued_barc COMPLETED!\n";
			}
		} elsif ($QUEUE{$queued_barc} == 0) { # store to kick off
			push @WAITING, $queued_barc;
		}
	}
	if ($running < $opt{'t'}) { # if not at full threading
		$ts = localtime(time);
#		print LOG "\t$ts\tThreads active: $running, adding...\n";
		for ($waitingID = 0; $waitingID < ($opt{'t'} - $running); $waitingID++) { # see how many short
			if (defined $QUEUE{$WAITING[$waitingID]}) { # if wiating is occupied
				system("$RealBin/sciMET_BSBolt2cellCalls.pl -s $opt{'O'} $WAITING[$waitingID] &"); # start the thread
				$QUEUE{$WAITING[$waitingID]} = 1; # log it as active
				print LOG "\t\t$WAITING[$waitingID] now running.\n";
			}
		}
	}
}
