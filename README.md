# nematodes

## to run exint-master.pl program on 120 Nematode species. It will get all the annotated gene sequences (exons-introns) using the genome and GFF file. 
        use strict;
        use warnings;

        my $filename = '/home/roylab/swadha/phylogeny/nematodes_V2/genenome_names.txt';
        open(my $fh, '<:encoding(UTF-8)', $filename)
          or die "Could not open file '$filename' $!";


        while (my $row = <$fh>) {
          chomp $row;
          #print "$row\n";
          my $x = $row;

          if (index($x, 'WBPS11.genomic.fa') != -1)
          {
            my $y = substr($x, 0, index($x, 'WBPS11.genomic.fa'));
            #print "$y\n"


            my $gff_ext = "WBPS11.annotations.gff3";
            my $genome_ext = "WBPS11.genomic.fa";



            my $gff_file = '/home/roylab/swadha/phylogeny/nematodes_V2/GFF_files/' . $y . $gff_ext;
            my $genome_file = '/home/roylab/swadha/phylogeny/nematodes_V2/genomes/' . $y . $genome_ext;


             print "$gff_file\n";
             print "$genome_file\n";

             system ("exint-master.pl $gff_file $genome_file > /home/roylab/swadha/phylogeny/nematodes_V2/step1_intron_exons/$y.exons-introns");
          }
          else
          {
            my $y = substr($x, 0, index($x, 'scaffolds.fa'));
            #print "$y\n"


            my $gff_ext = "gff3";
            my $genome_ext = "scaffolds.fa";



            my $gff_file = '/home/roylab/swadha/phylogeny/nematodes_V2/GFF_files/' . $y . $gff_ext;
            my $genome_file = '/home/roylab/swadha/phylogeny/nematodes_V2/genomes/' . $y . $genome_ext;


             print "$gff_file\n";
             print "$genome_file\n";

             system ("exint-master.pl $gff_file $genome_file > /home/roylab/swadha/phylogeny/nematodes_V2/step1_intron_exons/$y.exons-introns");

          }
        }

## exint-master.pl below
                        #!/usr/bin/perl

                        use English;
                        $maxline = 500;
                        open (IN, $ARGV[0]);



                        #  Get the CDs files, store information in %coords, %strand, and store the names of genes
                        #  on a given chromosome in %genes.
                        while (<IN>) {
                            ($chromosome,$lcoord,$rcoord,$strand) = (split)[0,3,4,6];
                            @a = split;
                            ($a[2] eq "CDS") || next;
                            #  Find the name
                            if (
                                /Parent=([^\s;]+)/ || 
                                /transcript_id "(\S+)"/ || 
                                /transcriptId (\S+)/ ||
                                /proteinId (\S+);/ || 
                                /protein_id="(\S+)";/ || 
                                /name "(\S+)"/i ||
                                /ID=([^\s;]+)/ || 
                                /mRNA\s+(\S+?)[\s\;]/ ||
                                /seq_id\s+(\S+)/
                                ) {
                                $gene = $1;
                            }
                            else {next;
                                #die "$ARGV[0] doesn't fit the format\n\n";
                            }
                            #  Add the current exon's coordinates
                            $coords{$gene} .= "$lcoord..$rcoord,";
                            #  Store the strand of the gene
                            $strand{$gene} = $strand;
                            #  Store the name of the gene in this array that is indexed on the chromosome
                            ($genes{$chromosome} =~ / $gene /) ||
                                ($genes{$chromosome} .= " $gene ");
                        }

                        #  If it didn't work for CDS's, it must be a file with 'exon' instead of 'CDS'.
                        #  So open the file up again and redo it that way.
                        unless (scalar keys %genes) {
                            #  Do it for the CDS-like file...
                            open (IN, $ARGV[0]);
                            while (<IN>) {
                                @a = split;
                                ($a[2] eq "exon") || next;
                                if (
                                    /CDS (\S+)/ || 
                                    /Parent=([^\s;]+)/ || 
                                    /transcript_id "(\S+)"/ || 
                                    /transcriptId (\S+)/ ||
                                    /ID=([^\s;]+)/ || 
                                    /mRNA (\S+);/) {
                                    $gene = $1;
                                }
                                else {
                                    die "$ARGV[0] doesn't fit the format\n\n";
                                }


                                $coords{$gene} .= "$a[3]..$a[4],";
                                $strand{$gene} = $a[6];
                                ($genes{$a[0]} =~ / $gene /) ||
                                    ($genes{$a[0]} .= " $gene ");
                            }
                        }

                        #  Here we go through the fasta of the chromosomes    
                        open (IN, $ARGV[1]);
                        $/ = ">";
                        while (<IN>) {
                            #  Get the name of the chromosome
                            ($chr) = /(\S+)/;
                            /\n/;
                            #  Get the sequence
                            @S= $POSTMATCH =~ /[A-Z]{1,$maxline}/ig;
                            $line = length $S[0];

                            #  For each gene on the chromosome...
                            foreach $gene ($genes{$chr} =~ /\S+/g) {
                                #  Get the coordinates
                                $coords = join (",", sort {$a<=>$b} $coords{$gene} =~ /\d+\.\.\d+/g);
                                #  Get the left- and right-most coordinates
                                ($l,$r) = ($coords =~ /\d+/g)[0,-1];
                                #  Get the sequence of that region
                                $seq = join ("", @S[($l-1)/$line..($r-1)/$line+1]);
                                $offset = ($l-1)%$line;
                                $seq =~ tr/A-Z/a-z/;
                                #  Change it all to lowercase
                                @G = ($seq =~ /./g)[$offset..$r-$l+$offset];
                                $coords=~ s/\d+/$&-$l/eg;
                                #  For each exon change it all to uppercase.
                                foreach $exon ($coords =~ /\d+\.\.\d+/g) {
                                    ($l,$r) = $exon =~ /\d+/g;
                                    $e = join ("", @G[$l..$r]);
                                    $e =~ tr/a-z/A-Z/;
                                    @G[$l..$r] = $e =~ /./g;
                                }
                                if ($strand{$gene} eq "-") {
                                    $seq = join ("", reverse @G);
                                    $seq =~ tr/ACGTacgt/TGCAtgca/;
                                }
                                else {$seq = join ("", @G);}

                                #  Print it out!
                                print ">$gene $strand{$gene}$chr:$coords{$gene}\n$seq\n"

                            }

                        }

## To run ei-to-exint.pl on 130 nematode species. It will delete the intron sequnces and write intron positions in the header.

        while IFS= read -r line

                   do
                        ei-to-exint.pl $line > /home/roylab/swadha/phylogeny/nematodes_V2/step2_ei-to-exint/$line

                   done < filenames.txt

## ei-to-exint.pl code below

                        #!/usr/bin/perl
                        use English;


                        $/ = ">";
                        while (<>) {
                            chomp;
                            ($name) = /(\S+)/;
                            /\n/;
                            $seq = join ("", $POSTMATCH =~ /\S+/g);

                            $ints = "";
                            #$seq =~ s/\s//g;
                            #$seq =~ s/\/.+?\///g;
                            #$seq =~ s/\\.+?\\//g;
                            $pos = 0;
                            while ($seq =~ /[A-Z]+/g) {
                                $pos += length $&;
                                $ints .= sprintf "%d.%d ", 1+$pos/3, $pos%3;
                            }
                            $ints =~ s/\d+\.\d $//;
                            print ">$name $ints\n";
                            print join ("", $seq =~ /[A-Z]+/g, "\n");
                        }





## translate.pl to translate the gene sequences created in the previous step. The translated sequences will have intron positions in thier header.


        while IFS= read -r line

                   do
                        translate.pl $line > /home/roylab/swadha/phylogeny/nematodes_V2/step3_translated_intron_exons/$line

                   done < filenames.txt
## Translate.pl below
                        #!/usr/bin/perl
                        %t = (
                            'TCA','S',
                            'TCC','S',
                            'TCG','S',
                            'TCT','S',
                            'TTC','F',
                            'TTT','F',
                            'TTA','L',
                            'TTG','L',
                            'TAC','Y',
                            'TAT','Y',
                            'TAA','_',
                            'TAG','_',
                            'TGC','C',
                            'TGT','C',
                            'TGA','_',
                            'TGG','W',
                            'CTA','L',
                            'CTC','L',
                            'CTG','L',
                            'CTT','L',
                            'CCA','P',
                            'CCC','P',
                            'CCG','P',
                            'CCT','P',
                            'CAC','H',
                            'CAT','H',
                            'CAA','Q',
                            'CAG','Q',
                            'CGA','R',
                            'CGC','R',
                            'CGG','R',
                            'CGT','R',
                            'ATA','I',
                            'ATC','I',
                            'ATT','I',
                            'ATG','M',
                            'ACA','T',
                            'ACC','T',
                            'ACG','T',
                            'ACT','T',
                            'AAC','N',
                            'AAT','N',
                            'AAA','K',
                            'AAG','K',
                            'AGC','S',
                            'AGT','S',
                            'AGA','R',
                            'AGG','R',
                            'GTA','V',
                            'GTC','V',
                            'GTG','V',
                            'GTT','V',
                            'GCA','A',
                            'GCC','A',
                            'GCG','A',
                            'GCT','A',
                            'GAC','D',
                            'GAT','D',
                            'GAA','E',
                            'GAG','E',
                            'GGA','G',
                            'GGC','G',
                            'GGG','G',
                            'GGT','G',
                            );
                        while (<>) {
                            chomp;
                            /\S/ || next;

                            ($seq) = /(\S+)/;
                            $head = (split)[1];
                            $seq =~ tr/Uu/Tt/;
                            $seq =~ tr/a-z/A-Z/;
                            $seq =~ s/[^A-Z]//g;

                            $pro = "";
                            while ($seq =~ /(.{1,3})/g) {
                                $pro .= trans($&);
                            }

                            $pro =~ /\S/ &&
                                print ">$head\n$pro\n";
                        }


                        sub trans () {
                            if (exists $t{$&}) {$t{$&}}
                            else {"X"}
                        }



## To make database of the translated sequences in the previous step

        while IFS= read -r line
                   do
                         mkdir /home/roylab/swadha/phylogeny/nematodes_V2/step4.2_database_for_blast_2/$line
                         makeblastdb -in $line -dbtype prot -out /home/roylab/swadha/phylogeny/nematodes_V2/step4.2_database_for_blast_2/$line/$line
                        

                   done < filenames.txt
## to perform blastp with c.elegans's histone H2As as querry

        while IFS= read -r line
                   do         
                        cd /home/roylab/swadha/phylogeny/nematodes_V2/step4.2_database_for_blast_2/$line
                        blastp -query /home/roylab/swadha/Diploscapter/NEMATODE/H4.fa -db /home/roylab/swadha/phylogeny/nematodes_V2/step4.2_database_for_blast_2/$line/$line -out /home/roylab/swadha/Diploscapter/NEMATODE/1_BLAST_results/H4/$line -outfmt "7 qacc sacc evalue qstart qend sstart send qlen slen length" -evalue 1e-10

                   done < /home/roylab/swadha/phylogeny/nematodes_V2/step4.2_database_for_blast_2/filename.txt

## Apply some filters in the blastoutput
        while IFS= read -r line

                   do 
                         cat /home/roylab/swadha/phylogeny/nematodes_V2/step4.3_blast_results/$line | grep "^[^#;]"  > /home/roylab/swadha/phylogeny/nematodes_V2/step4.4_filteredBlastout/$line.txt     
           

                   done < filenames.txt

## to fetch blast hit sequnces

        import pandas as pd
        arry = []
        filenames=[]
        with open ("/home/roylab/swadha/phylogeny/nematodes_V2/step4.4_filteredBlastout/all_in_one.txt", "r") as f:
            #with open ("fetched_sequences.txt", "w") as g:
                for line in f:
                    arry = line.split("\t")
                    #print (arry[1])
                    l = len(arry)
                    #print(l)
                    with open (arry[0],"r") as spFile:
                        #allLines = splLines in spFile
                        takeNextLine = 0
                        for spLines in spFile:
                            if (takeNextLine == 1):
                                print(arry[2] + "\t" + arry[1] + "\t" + arry[0])
                                print(spLines)
                                takeNextLine = 0

                            if (arry[2] in spLines):
                                takeNextLine = 1
                                #print(spLines)


## Run HOMOLOGOUS_INTRON_FILTER.PL (pasted below) on the fetched Blast hit sequnces
                        #!/usr/bin/perl
                        $window = 0;

                        use English;
                        #  Note: for very large introns, this program causes large amounts of intron sequence to be stored
                        #  in memory, which may lead to memory problems.  The current solution to this is to 

                        %translate = ('TCA','S','TCC','S','TCG','S','TCT','S','TTC','F','TTT','F','TTA','L','TTG','L','TAC','Y','TAT','Y','TAA','X','TAG','X','TGC','C','TGT','C','TGA','X','TGG','W','CTA','L','CTC','L','CTG','L','CTT','L','CCA','P','CCC','P','CCG','P','CCT','P','CAC','H','CAT','H','CAA','Q','CAG','Q','CGA','R','CGC','R','CGG','R','CGT','R','ATA','I','ATC','I','ATT','I','ATG','M','ACA','T','ACC','T','ACG','T','ACT','T','AAC','N','AAT','N','AAA','K','AAG','K','AGC','S','AGT','S','AGA','R','AGG','R','GTA','V','GTC','V','GTG','V','GTT','V','GCA','A','GCC','A','GCG','A','GCT','A','GAC','D','GAT','D','GAA','E','GAG','E','GGA','G','GGC','G','GGG','G','GGT','G');

                        ($#ARGV >= 1) || die "usage: HOMOLOGOUS_INTRON.PL blah.exons-introns blahblah.exons-introns.\n-window=length can be included anywhere on the commandline.\n\n";


                        $run = int (rand () * 1000000);

                        $fullstartover = $startover = 1;
                        @copy = @ARGV;

                        for ($i = $#ARGV; $i >= 0; $i--) {
                            if ($ARGV[$i] eq "-nostartover") {
                                $startover = 0;
                                $fullstartover = 0;
                                splice (@ARGV,$i,1);
                            }
                            elsif ($ARGV[$i] eq "-partstartover") {
                                $startover = 0;
                                splice (@ARGV,$i,1);
                            }
                            elsif ($ARGV[$i] eq "-firstversus") {
                                $firstversus = 1;
                                splice (@ARGV,$i,1);
                            }
                            elsif ($ARGV[$i] =~ /-orthologsfile=/) {
                                $orth_group_file_name = $POSTMATCH;
                                $startover = 0;
                                $fullstartover = 0;	
                                $use_provided_ortholog_list = 1;
                                splice (@ARGV,$i,1);
                            }
                            elsif ($ARGV[$i] =~ /-all/) {
                                $all = 1;
                                splice (@ARGV, $i, 1);
                            }
                            elsif ($ARGV[$i] =~ /-filter/) {
                                $filter = 1;
                                unless ($window) { $window = 10 }
                                unless ($cutoff) { $cutoff = 0.3 }
                                splice (@ARGV, $i, 1);
                            }
                            elsif ($ARGV[$i] =~ /-window=/) {
                                $window = $POSTMATCH;
                                splice (@ARGV,$i,1);
                            }
                            elsif ($ARGV[$i] =~ /-cutoff=/) {
                                $cutoff = $POSTMATCH;
                                splice (@ARGV,$i,1);
                            }
                            elsif ($ARGV[$i] =~ /-allowgaps/) {
                                $gaps = 1;
                                splice (@ARGV,$i,1);
                            }
                        }


                        if ($use_provided_ortholog_list) {
                            open (IN, $orth_group_file_name);
                            while (<IN>) {
                                foreach $g (/\S+/g) {
                                    $needed_genes{$g}++;
                                }
                            }
                        }



                        @sp = @ARGV;

                        #  Get data, make exint.pro files
                        foreach $sp (@sp) {
                            print STDERR "Here";
                            ($root{$sp}) = $sp =~ /(?:.+\/|)(\S*?)(\.|$)/;
                            open (IN, $sp);
                            unless ($use_provided_ortholog_list) {
                                open (OUT, ">$root{$sp}.exint.pro");
                                print STDERR "Making $root{$sp}.exint.pro\n";
                            }

                            $/ = ">";
                            <IN>;
                            while (<IN>) {
                                chomp;
                                ($name,$sequence) = /(\S+).*?\n(.+)/s;
                                (!$use_provided_ortholog_list) || $needed_genes{$name} || next;
                                ($sequence =~ /\S/) || next;

                                if ($sequence =~ /[\/\\]/) {
                                    $sequence =~ s/(\S+)\s+(\S+)\s+(\S+)/$1$2$3/g;
                                    @{$codons{$name}} = $sequence =~ /\S+/g;
                                }
                                else {
                                    $sequence = join ("", $sequence =~ /\S+/g);
                                    @{$codons{$name}} = $sequence =~ /[a-z]*[A-Z][a-z]*[A-Z][a-z]*[A-Z]/g;
                                }



                                $sequence =~ s/[[\/\\].+?[\/\\]//g;
                                $sequence =~ s/ //g;
                                while ($sequence =~ s/[^A-Z]+//) {
                                    $intron_positions{$name} .= sprintf "%d.%d ", $-[0]/3+1,$-[0]%3;
                                }
                                $sequence =~ s/[A-Z]{3}/$translate{$&}/g;
                                $pro{$name} = $sequence;
                                $use_provided_ortholog_list && next;
                                ($pro{$name} =~ /\S/) && 
                                    print OUT ">$name $intron_positions{$name}\n$pro{$name}\n";

                            }
                        }


                        unless ($use_provided_ortholog_list) {
                        print STDERR "Out of the exint.pro section...\n";
                        #  Do all-against-all blast searches in both directions
                        foreach $sp1 (@sp) {
                            system "formatdb -i $root{$sp1}.exint.pro";
                            foreach $sp2 (@sp) {
                                ($sp1 eq $sp2) && next;
                                print STDERR "In the blast section...\n";
                                (open (IN, "$root{$sp2}.v.$root{$sp1}") && !$fullstartover) && next;
                                print STDERR "blasting $root{$sp2}.exint.pro against $root{$sp1}.exint.pro\n";
                                system "blastall.pl -i $root{$sp2}.exint.pro -d $root{$sp1}.exint.pro -p blastp -e0.0000000001 -m8 > $root{$sp2}.v.$root{$sp1}";
                            }
                        }

                        #  Define ortholog pairs by best reciprocal hits

                        foreach $i (0..$#sp) {
                            foreach $j ($i+1..$#sp) {
                                print STDERR "Defining orthologs for $sp[$i] and $sp[$j].\n";
                                push (@file_list, $outfile = join (".v.", sort ($root{$sp[$i]},$root{$sp[$j]})) . ".orthologs");
                                unless ($use_provided_ortholog_list) {
                                    system "ORTHOLOG_PAIRS.PL $root{$sp[$i]}.v.$root{$sp[$j]} $root{$sp[$j]}.v.$root{$sp[$i]} > $outfile";
                                }
                                push (@file_list, (join (".v.", sort ($root{$sp[$i]},$root{$sp[$j]})) . ".orthologs"));    
                            }
                        }



                        close OUT;
                        }
                        ($orth_group_file_name =~ /\S/) || 
                            ($orth_group_file_name = join (".v.", sort values %root) . ".orthologs");


                        (@sp>2) &&
                            ($startover || !open (IN, $orth_group_file_name)) &&
                            system "ORTHOLOG_GROUPS.PL @file_list > $orth_group_file_name";



                        if ($all) {
                            $orth_group_file_name = "HI.ListOfAllGenes";
                            open (OUT, ">$orth_group_file_name");
                            foreach $g (keys %pro) {
                                print OUT "$g\t";
                            }
                            print OUT "\n";
                            close OUT;
                        }



                        #  Align and report ortholog sets
                        $/ = "\n";
                        open (ORTH, "$orth_group_file_name");
                        while (<ORTH>) {
                            $set_numb++;
                            (@genes) = /\S+/g;
                            open (OUT, ">clustin$run");
                            foreach $gene (@genes) {
                                print OUT ">$gene\n$pro{$gene}\n";
                            }
                            system "nohup clustalw clustin$run -output=gde -case=upper -outorder=input -quiet";
                            open (IN, "clustin$run.gde");
                            {
                                local $/ = "%";
                                $j = <IN>;
                                $g = 0;
                                while (<IN>) {
                                    chomp;
                                    $pos = 0;
                                    ($seq) = /\n(.+)/s;
                                    $name = $genes[$g];
                                    $seq =~ s/\s//g;

                                    @{$clust_aa{$name}} = $seq =~ /./g;
                                    $seq =~ s/-/--->/g;
                                    $seq =~ s/[A-Z]/$codons{$name}[$pos++].">"/eg;

                                    @{$clust_nt{$name}} = $seq =~ /[^>]+/g;
                                    $g++;
                                }
                            }
                            print ">$set_numb:" . join (",", @genes);
                            print "\nNumb:" . $#{$clust_aa{$genes[0]}} . "\n";


                            %id = ();
                            %conserved = ();
                            unless ($all) {
                                foreach $g1 (0..$#genes) {
                                    foreach $g2 ($g1+1..$#genes) {
                                        for $i (0..$#{$clust_aa{$genes[0]}}) {
                                            $t = $clust_aa{$genes[$g1]}[$i].$clust_aa{$genes[$g2]}[$i]; 
                                            if ($t eq "--") {$symbol = "G"}
                                            elsif ($t =~ /^-/) {$symbol = "F"}
                                            elsif ($t =~ /-/) {$symbol = "H"}
                                            elsif ($t =~ /(.)\1/) {$symbol = "1"}
                                            else {$symbol = 0}		    
                                            $id{$genes[$g1]}{$genes[$g2]}[$i] = 
                                                $id{$genes[$g2]}{$genes[$g1]}[$i] = $symbol;

                                        }
                                    }
                                }

                            }
                            for $i (0..$#{$clust_aa{$genes[0]}}) {
                                foreach $gene (@genes) {
                                    print "$clust_aa{$gene}[$i]\t";
                                }
                                foreach $gene (@genes) {
                                    print "$clust_nt{$gene}[$i]\t";
                                }
                                if (!$all && $window) {
                                    $no = 0;
                                    foreach $g1 (0..$#genes) {
                                        $firstversus && $g1 && last;
                                        $no && last;
                                        foreach $g2 ($g1+1..$#genes) {
                                            $t = join ("", reverse @{$id{$genes[$g1]}{$genes[$g2]}}[0..$i-1]);
                                            $t =~ /(.*?[01].*?){0,$window}[FGH]*/;
                                            $left = reverse $&;
                                            $pos = $id{$genes[$g1]}{$genes[$g2]}[$i];
                                            $t = join ("", @{$id{$genes[$g1]}{$genes[$g2]}}[$i+1..$#{$clust_aa{$genes[0]}}]);
                                            $t =~ /(.*?[01]){0,$window}[FGH]*/;
                                            $right = $&;

                                            if ($filter) {
                                                unless ($gaps || !(($left.$right) =~ /[FGH]/)) { $no++; last; }
                                                unless (($left =~ tr/01/01/) && ($right =~ tr/01/01/)) { $no++; last; }
                                                unless ($cutoff <= (($left =~ tr/1/1/) / ($left =~ tr/01/01/))) { $no++; last; }
                                                unless ($cutoff <= (($right =~ tr/1/1/) / ($right =~ tr/01/01/))) { $no++; last; }
                                            }
                                            else {
                                                print "$left($pos)$right:";
                                            }


                                        }
                                    }
                                    if ($filter) {
                                        print (("Y","N")[$no > 0]);
                                    }


                                }
                                print "\n";

                            }

                        }

## Fetch all the sequnces which has introns from the previous output
                i = 0
                with open ("/home/roylab/swadha/phylogeny/nematodes_V2/step6_groupings/codons.txt", "r") as fname:
                    for line in fname :
                        arry = line.split("\t")
                        i = i + 1
                        l = len(arry)
                        c1 = 1
                        c3 = 0
                        filename = "intron_pos_sequence_row_" + str(i) +  ".txt"
                        with open (filename, "w") as g:
                            for x in range (l):
                                columnNo = x
                                if any(c.islower() for c in (arry[x])):
                                    with open ("/home/roylab/swadha/phylogeny/nematodes_V2/step6_groupings/codons.txt", "r") as fname2:
                                        for line2 in fname2 :
                                            arry2 = line2.split("\t")
                                            g.writelines(arry2[x])
                                    g.writelines("\n")
                                    fname2.close()
                        g.close()
                fname.close()


## to fetch intronless sequneces
                import sys
                from collections import defaultdict

                datafile = sys.argv[1]

                columnwise = defaultdict(list)

                has_intron = set()

                with open(datafile) as data:
                    headers = data.readline().strip().split('\t')
                    for row in data:
                        columns = row.strip().split('\t')
                        for index, c in enumerate(columns):
                            if index in has_intron:
                                continue
                            name = headers[index]
                            if not (c.isupper() or c == '---'):  # is intron
                                has_intron.add(index)
                                if name in columnwise:
                                    del columnwise[name]
                                continue
                            # otherwise, no intron so add to column list
                            columnwise[name].append(c)

                for name, column_list in sorted(columnwise.items()):
                    column_list.insert(0, name)
                    print('\t'.join(column_list))


