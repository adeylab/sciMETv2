#!/usr/bin/perl

$die = "

sciMET_sortChroms.pl [chrom folder] (run)

Wrapper to sort the files in a chroms folder.
Will kick each off simultaneously using as many
threads as there are files.

Wrapper for linx sort command.
Sorts with -k1,1 -k2,2n, gzips and then deletes
original file, replacing it with the sorted one.

Prints the command to copy/paste.

If 'run' is specified as a second argumnet, will
run the command itself.

Once complete, remove the non gzipped bed files.

";

if (!defined $ARGV[0]) {die $die};

$command = "for f in $ARGV[0]/*.bed; do cat \$f | sort -k2,2n | gzip > \$f.gz & done";

print STDERR "Command:\n\t$command\n";

if (!defined $ARGV[1] || $ARGV[1] ne "run") {
	print STDERR "Copy/paste then run command.\n\n";
} elsif ($ARGV[1] eq "run") {
	print STDERR "Running!\n";
	system($command);
}
exit;