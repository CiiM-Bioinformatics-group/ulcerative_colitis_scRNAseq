my %hc;
my %sc;
my %use;

while(<>){
	chomp;
	next if (/^Gene/);
	my @it = split(/\t/);
	my @a = split(/\./,$it[4]);
	my $gene = $it[0];
	my $fc = $it[2];
	my $cell = $a[1];
	if($a[0] eq "UCvHC"){
		my $line = "$gene\t$cell";
		$hc{$line} = "$fc";
		$use{$line} = 1;
	}elsif($a[0] eq "UCvSC"){
		my $line = "$gene\t$cell";
		$sc{$line} = "$fc";
		$use{$line} = 1;
	}
}

print "Gene\tCell\tHC\tSC\n";
foreach my $key (keys %use){
	$hc{$key} = 0 unless($hc{$key});
	$sc{$key} = 0 unless($sc{$key});
	print "$key\t$hc{$key}\t$sc{$key}\n";
}



