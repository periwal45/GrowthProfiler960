#Program to parse stacker-reader output files in .txt format and generate tab separated files by plate.
#!/usr/bin/perl
use strict;
use warnings;

opendir(DIR, "/Users/periwal/ShikiFactory/WP3/GrowthProfiler/OD") || die $!;
my @dir = readdir(DIR);

foreach my $f1(@dir){

	next unless $f1 =~ /(^SF100.*)/;	#to read all bug folders starting with NT

	my $name1 = $1;

	opendir(NT, "/Users/periwal/ShikiFactory/WP3/GrowthProfiler/OD/$f1") || die $!;
	my @contents = readdir(NT);

	foreach my $f2(@contents){

		next unless $f2 =~ /(^all_bugs.*)/;	#to read all Replicate folders starting with Replicate

		my $name2 = $1;

		opendir(REP, "/Users/periwal/ShikiFactory/WP3/GrowthProfiler/OD/$f1/$f2") || die $!;
		my @rep = readdir(REP);

		foreach my $f3(@rep){

		next unless $f3 =~ /(^MTP.*)\.csv/;

		my $name3 = $1;

		#$name3 = lc($name3);

		my $file = "/Users/periwal/ShikiFactory/WP3/GrowthProfiler/OD/$f1/$f2/$f3";
		open my $fh, $file or die $!;
	
		my $outfile = "$name1"."_"."$name2"."_"."$name3".".tab";
		open(OUT, ">/Users/periwal/ShikiFactory/WP3/GrowthProfiler/OD/$f1/$outfile");

		while(my $line = <$fh>){

			#chomp $line;

			if ($line =~ /^Time \(min\)/) {
        		
        		$line =~ s/Time \(min\)/Time/g;
        		$line =~ s/,,.*//g;
        		$line =~ s/,/\t/g;

        		print OUT "$line";

        		}elsif($line =~ /^\d+,/){

        			$line =~ s/,,.*//g;
        			$line =~ s/,/\t/g;
        		
        		print OUT "$line";

        		}

        		

        	
   			 } 
		}	
	}
}
	

close DIR;
close NT;
close REP;
close OUT;	
