#!/usr/bin/perl
#
#############################################################################
# Program       : clips                                                     #
# Author        : Luis Mendoza                                              #
# Date          : 07.07.18                                                  #
# SVN Info      : $Id: clips.pl 8811 2023-01-10 22:55:10Z real_procopio $
#                                                                           #
# Generate reverse index of n-AA keys, for use in mapping observed          #
#  peptide sequences to all proteins                                        #
# Copyright (C) 2018-2023 Luis Mendoza                                      #
#                                                                           #
# This library is free software; you can redistribute it and/or             #
# modify it under the terms of the GNU Lesser General Public                #
# License as published by the Free Software Foundation; either              #
# version 2.1 of the License, or (at your option) any later version.        #
#                                                                           #
# This library is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         #
# General Public License for more details.                                  #
#                                                                           #
# You should have received a copy of the GNU Lesser General Public          #
# License along with this library; if not, write to the Free Software       #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
#                                                                           #
# Institute for Systems Biology                                             #
# 401 Terry Avenue North                                                    #
# Seattle, WA  98109  USA                                                   #
# lmendoza@isb                                                              #
#                                                                           #
#############################################################################
#
use strict;
use warnings;
use Getopt::Std;
use Tie::File;
use FindBin qw($Bin);
use lib "$Bin/../lib/perl";

my $VersionInfo = "stand-alone 1.5";
my $rc = eval {
  require tpplib_perl; # exported TPP lib function points
  tpplib_perl->import();
  1;
};
if ($rc) { 
  $VersionInfo = tpplib_perl::getTPPVersionInfo();
}

$|++; # autoflush

my $USAGE=<<"EOU";
--------------------------------------------------------------------------
Program:    $0 ($VersionInfo)
Purpose:    Generates protein sequence index file by use of segments, for
            use in mapping observed peptide sequences to all proteins.
            Works with any protein file in FASTA format, including PEFF.

Usage:      $0 [options] <fasta_file>
Generates:  <fasta_file>.pep.idx
Options:
            -s <length>  segment size, in number of amino acids [default=5]
            -V           do not use PEFF variants [default:use them]
            -f           force overwriting of index file, if exists
            -I           do NOT convert I->L
            -X <string>  exclude protein entries that start with <string>
                         (e.g. DECOY_)
            -A           do not generate all possible keys in index
                         (not recommended for large files)
For Developers:
            -D           print debug information
--------------------------------------------------------------------------
EOU

my %options;
getopts('s:VfIADX:', \%options);

my $DEBUG      = $options{'D'} || '';
my $inpeff     = shift || die "Error: Please provide an input fasta file.\n".$USAGE;
my $outfile    = "$inpeff.pep.idx";
my $pepkeysize = $options{'s'} || 5;
my $excludestr = $options{'X'} || '---@@@---***---@@@---';
my $usepeff    = $options{'V'} ? 0:1;

my $aa_list    = 'ACDEFGHIKLMNPQRSTVWY';
$aa_list =~ s/I// unless $options{'I'};

print
    "Creating index for $inpeff; using keys of length $pepkeysize.\n",
    "(". ($usepeff?'':'NOT ') ."using PEFF variants)",
    "(". ($options{'I'}?'NOT ':'') ."mapping I->L)",
    "(". ($options{'A'}?'NOT ':'') ."creating full index keys)",
    $options{'X'} ? "(excluding '$excludestr' entries)" : '',
    $DEBUG ? '(debug mode)' : '',
    "\n";

!((-e $outfile) && !$options{'f'}) || die "\n\nIndex file $outfile already exists!  Delete or rename the file, or use the -f option to overwrite.  Exiting";

# test that we can do this so we can give error early on...
open(INDX, ">$outfile") || die "\n\nError: cannot open $outfile. Exiting";
close (INDX);

my %index;
my %segidx;
my %prots;
my %protseqs;

&processFile('enter');
&getAliases();
&processFile('prots');
&writeIndex();

exit(0);


####################################################################################
sub processFile {
####################################################################################
    my $do = shift || 'prots';

    print "Calculating protein lengths..." unless $do eq 'prots';

    my $protseq = '';
    my $prothdr = '';
    my $i = 0;

    open(PEFF, $inpeff) || die "\n\nError: cannot open $inpeff. Exiting";
    while(<PEFF>) {
	next if (/^#/);
	s/[\r\n]//g; # cheap chop

	if (/^>/) {
	    $i++;
	    if ($do eq 'prots') {
		# progress status...
		print (" $i\t") unless ($i % 100);
		print ("\n") unless ($i % 1000);

		&processProt($prothdr,$protseq);
	    }
	    else {
		&addProt22($prothdr, $protseq);
#		&addProt($prothdr, length $protseq);
	    }
	    $prothdr = $_;
	    $protseq = '';
	}
	else {
	    s/I/L/g unless $options{'I'};
	    $protseq .= $_;
	}

    }
    # get the last one
    if ($do eq 'prots') {
	&processProt($prothdr,$protseq);
	print "\n";
    }
    else {
	&addProt22($prothdr, $protseq);
#	&addProt($prothdr, length $protseq);
	print "ok.\n";
	print "Found " , scalar keys %protseqs , " unique protein sequences ($i total)\n";
	%protseqs = ();  # free memory; re-use to keep track of processed sequences
    }
    close (PEFF);

}


####################################################################################
sub processProt {
####################################################################################
    my $header = shift || return;
    $header =~ /^>(\S*)\s*(.*)/;
    my $acc = $1;
    my $annot = $2;

    my $protseq = shift;# || die "\n\nno sequence found for $acc.  Exiting...";

#    if (length($protseq) < $pepkeysize) {
#	print "\nWarning: protein $acc is too short. Skipping...\n";
#	return;
#    }

    my $alias = defined $prots{$acc} ? $prots{$acc} : $prots{"___alias___$acc"};

    return unless defined $alias;
    return if $protseqs{$alias};
    $protseqs{$alias} = 1;  # (pre-)mark as processed

    # process variants
    if ($annot =~ /VariantSimple=(\S*)/ && $usepeff) {
	my $v = $1;
	$v =~ s/I/L/g unless $options{'I'};

	my @protvar = ('', split("", $protseq));  # make it 1-based for speed and ease of variant processing

	# capture PEFF variants
	while ($v =~ m/\((\d+)\|(.)[\|\)]/g) {
	    $protvar[$1] .= $2;
	}


	my $new = 0;


	if ($new) {
	    my @variants = @protvar[1..$pepkeysize];
	    my %ret = ();
	    my %ppp;

	    &assemble_variants_NEW('', $alias, 1, \@variants, \%ret);

	    my $pos = 1;

	    for (my $i = $pepkeysize+1; $i <= length($protseq); $i++) {
		$pos++;

		%ppp = %ret;
		%ret = ();

		for my $k (keys %ppp) {
		    for my $aa (split("", $protvar[$i])) {
			next if $aa eq '*';
			$index{"$k$aa"} .= ":$alias,$pos";    # format:  PROT,POS

#			print "==>",substr("$k$aa",-($pepkeysize-1)),"<==\n" if ($i-$pepkeysize+1) == 2;

			$ret{substr("$k$aa",-($pepkeysize-1))} = 1;
		    }
		}
	    }

	}


	else {
	    for (my $i = 1; $i <= length($protseq)-$pepkeysize+1; $i++) {
		my @variants = @protvar[$i..$i+$pepkeysize];
		&assemble_variants('', $alias, $i, \@variants);
	    }
	}


    }
    else {
	for (my $i = 0; $i <= length($protseq) - $pepkeysize; $i++) {
	    my $key = substr $protseq, $i, $pepkeysize;
	    $index{$key} .= ":$alias,".($i+1);    # format:  PROT,POS
	    ## $index{$key} .= ":".($i+1).".$alias"; # format:  POS.PROT (allows numerical treatment)
	}
    }

}




####################################################################################
sub assemble_variants_NEW {    # recurse me...recurse me, my friend...
####################################################################################
    my $str = shift;
    my $prt = shift;
    my $pos = shift;
    my $varaas = shift;
    my $prestr = shift;

    for my $aa (split("", @$varaas[length($str)])) {
	next if $aa eq '*';
	if (1+length($str) == $pepkeysize) {
	    $index{"$str$aa"} .= ":$prt,$pos";    # format:  PROT,POS    #### speed this up with arrays??



	    ###$index{"$str$aa"} .= ":$prt,".($pos+1);    # format:  PROT,POS    #### speed this up with arrays??
	    ## $index{"$str$aa"} .= ":".($pos+1).".$prt"; # format:  POS.PROT (allows numerical treatment)
	    ## $ret{"$str$aa"} = 1;

#	    print "==>",substr("$str$aa",-($pepkeysize-1)),"<==\n";

	    $$prestr{substr("$str$aa",-($pepkeysize-1))} = 1;

#	    print "==$str$aa\n";



	}
	else {
	    &assemble_variants_NEW("$str$aa", $prt, $pos, $varaas, $prestr);
	}
    }
}







####################################################################################
sub assemble_variants {    # recurse me...recurse me, my friend...
####################################################################################
    my $str = shift;
    my $prt = shift;
    my $pos = shift;
    my $varaas = shift;

    for my $aa (split("", @$varaas[length($str)])) {
	next if $aa eq '*';
	if (1+length($str) == $pepkeysize) {
	    $index{"$str$aa"} .= ":$prt,$pos";    # format:  PROT,POS    #### speed this up with arrays??

### OLD	    $index{"$str$aa"} .= ":$prt,".($pos+1);    # format:  PROT,POS    #### speed this up with arrays??


	    ## $index{"$str$aa"} .= ":".($pos+1).".$prt"; # format:  POS.PROT (allows numerical treatment)
	    #print "$str$aa\n";
	}
	else {
	    &assemble_variants("$str$aa", $prt, $pos, $varaas);
	}
    }
}


####################################################################################
sub writeIndex {
####################################################################################
    # remove non-prots for statistics
    my @nonprots = grep {/___alias___/} keys %prots;
    delete @prots{@nonprots};

# report that we ignore * ?

    open(INDX, ">$outfile") || die "\n\nError: cannot open $outfile. Weird...";
    binmode INDX;

    print "Writing file...";

    my $header =
	"# Index generated by $0\n".
	'# Date=' . scalar(localtime) . "\n".
	"# OriginalFile=$inpeff\n".
        '# AASub=' . ($options{'I'} ? 'NONE': 'I->L') . "\n".
	"# KeyLength=$pepkeysize\n".
        '# PEFFVariants=' . ($usepeff ? 'VariantSimple':'NONE') . "\n".
	'# NumMappedProteinSequences=' . scalar(keys %prots) . "\n".
	'# NumSegments=' . scalar(keys %index) . "\n".
	"# BeginSegmentOffsets\n";

    print INDX $header;

    # placeholder for segments offset index
    my $zero = sprintf("%010d", 0);  # should be good to 10Gb!
    for my $a1 (split("", $aa_list)) {
	for my $a2 (split("", $aa_list)) {
#	    for my $a3 (split("", $aa_list)) {
	    my $a3 = '';
		print INDX "$a1$a2${a3}::$zero\n";
#	    }
	}
    }

    print INDX "# BeginProteins\n";
    # this is a sort of convenience
    for(sort keys %prots) {
	print INDX "$prots{$_}::$_\n" unless /^___alias___/;;
    }
    print INDX "# BeginIndex\n";

    # write all possible aa-combos in alpha order to index; set to empty if not found in hash (unless -A)
    &writeSegments('');
    close (INDX);

    # write segments offset index
    my $sis = () = $header =~ /\n/g;
    print "...seg-index...";

    # tie me...
    tie my @indexfile, "Tie::File", $outfile, recsep => "\n", memory => 0;
    for(sort keys %segidx) {
	$indexfile[$sis++] = sprintf("${_}::%010d", $segidx{$_});  # should be good to 10Gb!
    }
    untie @indexfile;

    print "finished. Index is at $outfile\n";
}

####################################################################################
sub writeSegments {    # recurse me...recurse me, again
####################################################################################
    my $kstr = shift;
    print $kstr if length($kstr) == 1;

    for my $aa (split("", $aa_list)) {
	if (length($kstr) == 1) {  # 2 for 3-AA
	    my $offset = tell INDX;
	    $segidx{"$kstr$aa"} = sprintf("%010d", $offset);
	}

	if (1+length($kstr) == $pepkeysize) {
	    if ($index{"$kstr$aa"}) {
		print INDX "$kstr$aa".$index{"$kstr$aa"}."\n";
	    }
	    elsif (!$options{'A'}) {
		print INDX "$kstr$aa:\n";
	    }
	}
	else {
	    &writeSegments("$kstr$aa");
	}
    }
}


####################################################################################
sub getAliases {
####################################################################################
    my @pkeys = reverse sort {(split /___/,$a)[0]<=>(split /___/,$b)[0]} keys %prots;
    %prots = ();

    my @alphabet = (0..9,'a'..'z','A'..'Z'); # cheap base-62 compression of aliases

#    print "Found ". scalar @pkeys ." protein sequences.\n";
    for my $i (0..$#pkeys) {
	next if $pkeys[$i] =~ /^00___uniq___/;

	$pkeys[$i] =~ /.*___(.*)/;

	my $plist = $1;

	my $rem = $i % scalar @alphabet;
	my $int = int($i / scalar @alphabet);

	my $b62 = "$alphabet[$rem]";

	while ($int > 0) {
	    $rem = $int % scalar @alphabet;
	    $b62 = "$alphabet[$rem]$b62";

	    $int = int($int / scalar @alphabet);
	}

	$prots{$plist} = $b62;

	print "$pkeys[$i]\t\t$prots{$plist}\t\t$plist\n" if $DEBUG;

	if ($plist =~ /\s/) {
	    for my $p (split ' ', $plist) {
		$prots{"___alias___$p"} = $b62;
	    }
	}
    }
}

####################################################################################
sub addProt {
####################################################################################
    my $header = shift || return;
    $header =~ /^>(\S*)\s*/;
    my $acc = $1;

    my $protlen = shift || die "\n\nno length found for $1.  Exiting...";

    $prots{"${protlen}___$acc"} = 1;
}
####################################################################################
sub addProt22 {
####################################################################################
    my $header = shift || return;
    $header =~ /^>(\S*)\s*/;
    my $acc = $1;

    if ($acc =~ /^\Q$excludestr\E/) {
	print "Info: protein $acc is excluded as per user input. Skipping...\n" if $DEBUG;
	return;
    }

    #    my $sequence = shift || die "\n\nno sequence found for $1.  Exiting...";

    my $sequence = shift;
    if (!$sequence) {
      print "Warning: no sequence found for protein $acc. Skipping...\n";
      return;
    }
    my $protlen = length $sequence;

    if ($protlen < $pepkeysize) {
	print "Warning: protein $acc is too short. Skipping...\n";
	return;
    }

    my $alreadyseen = 0;

    if ($prots{"00___uniq___$acc"}) {
	$alreadyseen = 1;
    }
    else {
	$prots{"00___uniq___$acc"} = 1;
    }


    if ($header =~ /VariantSimple=(\S*)/ && $usepeff) {
	$sequence = $1.$sequence;
    }

    if ($protseqs{$sequence}) {
	print "Found identical sequences!  $protseqs{$sequence} <--> $acc\n" if $DEBUG;

	if ($alreadyseen) {
	    for my $p (split ' ', $protseqs{$sequence}) {
		$alreadyseen*= ($acc ne $p);
	    }
	}

	delete $prots{"${protlen}___$protseqs{$sequence}"};

	$protseqs{$sequence} .= " $acc";
    }
    else {
	$protseqs{$sequence} = $acc;
    }

    print "Warning!  More than one entry for $acc with DIFFERENT sequences! Please check your database.\n" if $alreadyseen;

    $prots{"${protlen}___$protseqs{$sequence}"} = 1;
}
