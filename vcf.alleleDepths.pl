use strict;
use warnings;
use Getopt::Long;

GetOptions('p' => \(my $pass = ''), 'n=s' => \(my $normalGenotype = ''));
my ($vcfFile, $normalSample, @tumorSampleList) = @ARGV;
my @sampleList = ();
open(my $reader, $vcfFile);
while(my $line = <$reader>) {
	chomp($line);
	next if($line =~ s/^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t// && (@sampleList = split(/\t/, $line)));
	next if($line =~ /^#/ || scalar(@sampleList) == 0);
	my ($chromosome, $position, $id, $refBase, $altBase, $quality, $filter, $info, $format, @genotypeList) = split(/\t/, $line);
	next unless($pass || $filter eq 'PASS');
	my %sampleGenotypeHash = ();
	for(my $index = 0; $index < scalar(@sampleList); $index++) {
		@{$sampleGenotypeHash{$sampleList[$index]}}{split(/:/, $format)} = split(/:/, $genotypeList[$index]);
	}
	next if($normalGenotype ne '' && $normalGenotype ne $sampleGenotypeHash{$normalSample}->{'GT'});
	next if(grep {$_->{'GT'} eq './.'} @sampleGenotypeHash{$normalSample, @tumorSampleList});
	my @alleleDepthListList = map {[split(/,/, $_->{'AD'})]} @sampleGenotypeHash{$normalSample, @tumorSampleList};
	next if(grep {scalar(@$_) != 2} @alleleDepthListList);
	print join("\t", $chromosome, $position, $refBase, $altBase, map {@$_} @alleleDepthListList), "\n";
}
close($reader);
