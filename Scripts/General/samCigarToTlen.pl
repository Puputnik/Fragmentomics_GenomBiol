#!/usr/bin/env perl
use strict;


LINE: while (my $l = <STDIN>)
  {
    chomp $l;
    my @f=split(/\t/,$l);

    if ( ($f[0] !~ /^\@/) &&
         ($f[2] ne '*') )
      {
        my $len=$f[5];

        # Code from filippo
        my @cig = ( $f[5] =~ /(\d+\w)/g);
        my $len = 0;
        for my $c (@cig) {
          if ($c =~ /(\d+)[MX=DN]/g) {  # This takes the length of the reference, not the query (because insertions are usually wrong in nanopore)
            $len = $len + $1;
          }
        }

        $f[8]=$len;
      }
    print join("\t",@f)."\n";
  }
