while(<>){
	chomp;
	s/[\,\s]+/\n/g;
	print "$_\n";
}

