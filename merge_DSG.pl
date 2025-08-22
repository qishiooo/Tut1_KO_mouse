use strict;
use File::Basename;


open OUT,">DSG.txt" or die $!;
print OUT "Group\tType\tGene_name\tPvalue\tdPSI\tCoordinate\tWT\tPKO\tTKO\tTPDKO\tIsoform_IDs\n";
for my $file(glob "diffsplice_all/*/*_vs_*.*.txt"){
	my $name=basename($file,".txt");
	my($group, $type)=split /\./,$name;
	if($type eq 'es'){
		$type='AE';
	}elsif($type eq 'ir'){
		$type='IR';
	}elsif($type eq 'alt3'){
		$type='A3SS';
	}elsif($type eq 'alt5'){
		$type='A5SS';
	}
	open IN,"$file" or die $!;
	<IN>;
	while(<IN>){
		chomp;
		my @F=split /\t/;
		my $out=join "\t", $group, $type, @F;
		if( abs($F[2])>=10 ){
			print OUT "$out\n";
		}
	}
	close IN;
}
close OUT;



