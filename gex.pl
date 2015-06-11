
die "genes.txt, expression data.txt\n" unless (@ARGV);
my ($gf,$ef)=@ARGV;

my %GLIST; print STDERR read_genelist(\%GLIST,$gf), " genes found in $gf\n";
my $ln=0; $nfound=0;
open (FILE, "<",$ef) or die "$ef not read\n";
while (my $line=<FILE>) {
        if ($ln==0) { print $line; }
        else {
                my @vals=split "\t", $line;
                if ($GLIST{$vals[0]} ne "") { print $line; $nfound++; }
        }
        $ln++;
}
close FILE;
print STDERR $ln, " lines read $nfound were matched\n";


sub read_genelist {
        my ($data,$file)=@_;
        open (FILE, "<", $file) or return 0;
        while (my $line=<FILE>) { chomp $line; $line=~s/\s//g;
                $$data{$line}++;
        }
        close FILE;
        return scalar keys %{$data};

}

