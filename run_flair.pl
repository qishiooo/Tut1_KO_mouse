use strict;
use my_utils;
use Data::Dumper;
use File::Basename;

##  conda activate flair


my $ref="/data/Data/database/mmu/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa";
my $gtf="/data/Data/database/mmu/Mus_musculus.GRCm38.102.gtf";


my @samples=qw /WT-1 WT-2 WT-3 PKO-1 PKO-2 PKO-3 TKO-1 TKO-2 TKO-3 TPDKO-1 TPDKO-2 TPDKO-3/;
&my_threads(12, \&run, \@samples);
system "cat *.correct_all_corrected.bed  >correct/all.correct.bed";

my @chrs= (1..19, 'MT', 'X', 'Y');
my @reads= map { "../full/$_.full.fq.gz" } @samples;
&my_threads(3, \&collapse, \@chrs);

system "cat collapse/all.chr*.collapse.isoforms.fa >collapse/all.isoforms.fa";
system "flair quantify --threads 80 -r reads_manifest.tsv -i collapse/all.isoforms.fa --output quantify/all.quantify --temp_dir temp ";
system "cat all.chr*.collapse.isoforms.bed |sort -u  > collapse/all.isoforms.bed";
system "mark_producivity collapse/all.isoforms.bed $gtf $ref >productivity.bed";
system "mark_intron_retention productivity.bed mark_IR_isoforms.bed  mark_IR_introns.bed"; 

open IN,"group.txt" or die $!;
my @conA;
my @conB;
while(<IN>){
	chomp;
	my @F=split /\t/;
	push @conA, $F[0];
	push @conB, $F[1];
}
close IN;

open IN, "$gtf" or die $!;
my %genes;
while(<IN>){
	chomp;
	next if /^#/;
	my @F=split /\t/;
	next unless($F[2] eq 'gene');
	(my $id)=$F[-1]=~/gene_id "([^"]+)"/;
	(my $g)= $F[-1]=~/gene_name "([^"]+)"/;
	$genes{$id}= $g // '';
}
close IN;

print "@conA\n@conB\n";
&my_threads(5, \&diffsplice, \@conA, \@conB);
&diffsplice('TKO', 'WT');

sub run{
	my $s=shift @_;
	system "pychopper -r $s.report.pdf -u $s.unclassified.fq -w $s.rescued.fq $samples{$s}  $s.fq";
	system "flair align  -g $ref -r ../full/$s.full.fq.gz --output align/$s.aligned";
	my %ngs;
	for my $file(glob "/rawdata/*.R1.raw.fastq.gz"){
		my $sam=basename($file,".R1.raw.fastq.gz");
		(my $file2=$file)=~s/\.R1\./\.R2\./;
		(my $s)= $sam=~/\-(.+)/;
		$s=~s/_/-/;
		$ngs{$s}= ["$file", "$file2"];
	}
	system "/home/software/hisat2/hisat2/hisat2 -p 4 --dta -x $ref -1 $ngs{$s}->[0]  -2 $ngs{$s}->[1] -S hisat2/$s.sam";
	system "junctions_from_sam -s hisat2/$s.sam -n $s";
	system "flair correct -q align/$s.aligned.bed  -f $gtf -g $ref  --shortread  $s\_junctions.bed  --output correct/$s.correct";
}

sub collapse{
	my $chr=shift @_;
	print "chr:$chr\n";
	system "grep '^$chr\\s' correct/all.correct.bed >correct/all.chr$chr.bed";
	system "flair collapse --threads  30 -g $ref --gtf $gtf -q correct/all.chr$chr.bed -r @reads --stringent --check_splice --generate_map --annotation_reliant generate  --output collapse/all.chr$chr.collapse  --temp_dir temp";
}

sub diffsplice{
	my($s1,$s2)=@_;
	my $dir="diffsplice_all/$s1\_vs_$s2";
	system "flair diffSplice --threads 20  -i collapse/all.isoforms.bed -q quantify/all.quantify.counts.tsv --conditionA $s1 --conditionB $s2 --out_dir $dir";
        for my $type(qw /alt3 alt5 es ir/){
                open IN,"$dir/diffsplice.$type.events.quant.tsv" or die $!;
                open OUT,">$dir/diffsplice.$type.events.quant.txt" or die $!;
                my %counts;
                while(<IN>){
                        chomp;
                        my @F=split /\t/;
                        if($.==1){
                                print OUT "feature_id\tcoordinate\tWT\tPKO\tTKO\tTPDKO\tisoform_ids\n";
                        }else{
                                my $wt=$F[11]+$F[12]+$F[13];
                                my $tko=$F[5]+$F[6]+$F[7];
                                my $dko=$F[8]+$F[9]+$F[10];
                                my $pko=$F[2]+$F[3]+$F[4];
                                my $out=join "\t",$F[0],$F[1],$wt,$pko,$tko,$dko,$F[-1];
								my $type=$F[0]=~/^inc/ ? 'inclusion' : 'exclusion';
                                $counts{$F[1]}{$type}={
                                        wt => $wt,
                                        tko => $tko,
                                        dko => $dko,
                                        pko => $pko,
                                };
                                print OUT "$out\n";
                        }
                }
                close IN;
                close OUT;
                system "diffsplice_fishers_exact $dir/diffsplice.$type.events.quant.txt $s1 $s2 $dir/fishers.$type.$s1\_vs$s2.txt";
                my %stat;
                my %pvalue;
                open IN,"$dir/fishers.$type.$s1\_vs$s2.txt" or die $!;
                my $head=<IN>;
                chomp $head;
                while(<IN>){
                        chomp;
                        my @F=split /\t/, ;
                        next unless $F[-1]<=0.001;
                        (my $id)=$F[-2]=~/(ENSMUSG\d+)/;
                        my $g='';
                        if($id and $genes{$id}){
                                $g=$genes{$id};
                        }
                        $pvalue{$F[-1]}{$F[1]}++;
						my $type=$F[0]=~/^inc/ ? 'inclusion' : 'exclusion';
                        $stat{$F[1]}{$type}={
                                p => $F[-1],
                                gene => $g,
                                iso => $F[-2],
                        };
                }
                close IN;
                open OUT,">$dir/$s1\_vs_$s2.$type.txt" or die $!;
                print OUT "Gene_name\t$s1-$s2\_pval\tdPSI\tcoordinate\tWT\tPKO\tTKO\tTPDKO\tisoform_ids\n";
                for my $p(sort {$a<=>$b} keys %pvalue){
                        for my $pos(keys %{$pvalue{$p}}){
								my ($g,$iso,%psi,$n1,$n2,$n3,$n4,);
                                for my $in(qw/inclusion exclusion/){
                                    $g=$stat{$pos}{$in}{gene} if !$g;
                                    $iso.= $iso ? ";$stat{$pos}{$in}{iso}" : $stat{$pos}{$in}{iso};
									$n1.= $in eq 'inclusion' ? $counts{$pos}{$in}{wt} : ";$counts{$pos}{$in}{wt}";
									$n2.= $in eq 'inclusion' ? $counts{$pos}{$in}{pko} : ";$counts{$pos}{$in}{pko}";
									$n3.= $in eq 'inclusion' ? $counts{$pos}{$in}{tko} : ";$counts{$pos}{$in}{tko}";
									$n4.= $in eq 'inclusion' ? $counts{$pos}{$in}{dko} : ";$counts{$pos}{$in}{dko}";
									$psi{$in}={
										'WT' => $counts{$pos}{$in}{wt},
										'PKO' => $counts{$pos}{$in}{pko},
										'TKO' => $counts{$pos}{$in}{tko},
										'TPDKO' => $counts{$pos}{$in}{dko},
									};
                                }
								my $dpsi= 100*($psi{inclusion}{$s1}/($psi{inclusion}{$s1}+$psi{exclusion}{$s1})-$psi{inclusion}{$s2}/($psi{inclusion}{$s2}+$psi{exclusion}{$s2}));
                                my $out=join "\t", $g, $p, $dpsi, $pos, $n1,$n2,$n3,$n4, $iso;
                                print OUT "$out\n";
                        }
                }
                close OUT;
        }
}

	



