#! /usr/bin/perl -w
use strict;
use IO::Handle;
STDOUT->flush();
my @samples = qw(AR5361  BB4911  BB4913  CD5941 DN4461  EK3261  EP5591  ET2371  JC5733  JFH3411  JM1191  JZ1851 LK4881  MA5421  MA5731  MC3481  QW4591  SL2221  SM2701  TA2691  TF3471);
#my @samples = qw(JFH3411 JM1191 JZ1851 KP6951 LK4881 MA5421);
#my @samples = qw(DI6941 DN4461 EK3261 EP5591 ET2371 JC5733);
#my @samples = qw(MA5731 MC3481 Penn342 QW4591 SL2221 SM2701 TA2691 TF3471);
system("mkdir -p /nfsscratch/eshaeffer/workspace/dna/vcfs/deep_variant_bedFilter");
my $sample;
my @somatic_filenames_array;
my @filenames_array;
my %samples_hash =();
my $size_of_filenames_array;
my @my_ar;
my $c;
my $num;
my $max_size_of_filename_array = 0;
my $fi;
my $f;
foreach $sample (@samples) {
   chomp $sample;
   @filenames_array = ();
   @somatic_filenames_array = ();
   print "FILENAMES $sample\n";  
   @filenames_array = `./find_filenames_both.sh $sample`;
   print $sample."\n";
   print $sample." GERMLINE AND SOMATIC  FILENAMES: ".scalar(@filenames_array)."\n";
   print $filenames_array[0]."\n";
   $fi = 0;
   foreach $f (@filenames_array) {
      $samples_hash{$sample}{$fi} = $f;
      $fi++;
   }
   $size_of_filenames_array = scalar(@filenames_array);
   if ($size_of_filenames_array > $max_size_of_filename_array) {
      $max_size_of_filename_array = $size_of_filenames_array;
   }
}
my $filename = "";
my $num_filtered = 0;
my $num_passes = 0;
my $OUT;
my $PASSES;
my $file_prefix = "redo_";
foreach $sample (@samples) {
   chomp $sample;
   $file_prefix .= $sample."_";
}
open $OUT, ">", $file_prefix."nums_bed_filtered_samples.txt" or die "Cannot open nums_bed_filtered_saples.txt:$!\n";
open $PASSES, ">", $file_prefix."nums_bed_filtered_passed_samples.txt" or die "Cannot open nums_bed_filtered_saples.txt:$!\n";
foreach $sample (@samples) {
   chomp $sample;
   print $OUT $sample."\t";
   print $PASSES $sample."\t";
}
print $OUT "\n";
print $PASSES "\n";

my $ind;
for $ind (0 ...($max_size_of_filename_array - 1)) {
   print "SAMPLE F ___________\n";
   print "SAMPLE F index: ".$ind."\n";
   $filename = "";
   foreach $sample (@samples) {
      print "SAMPLE F $sample\n";
      $num_filtered = 0;
      print "sample: ".$sample."\n";
      #$filename = $samples_hash{$sample}->[$ind];
      $filename = $samples_hash{$sample}{$ind};
      print "SAMPLE F $filename\n";
      if ($filename =~ m/\w.*vcf/) {
          print "SAMPLE F $sample contains filename\n";
          chomp $filename;
          $num_filtered =  run_bed_filter_on_file($filename);
          $num_passes = $num_filtered;
          $num_filtered =~ s/^(.*)\t.*$/$1/;
          $num_passes =~ s/^.*\t(.*)$/$1/;

          print $OUT $num_filtered; 
          print $PASSES $num_passes;
          print "SAMPLE F: $sample ".$num_filtered."\n"; 
      } else {
          print $OUT "";
          print $PASSES "";
          print "";
      }
      #print "\t";
      print $OUT "\t";
      print $PASSES "\t";
   }
   print $OUT "\n";
   print $PASSES "\n";
   print "\n";
}


sub run_bed_filter_on_file {
   my $file = shift;
   chomp $file;
   my $count_pass = 0;
   my $outdir = "/nfsscratch/eshaeffer/workspace/dna/vcfs/bedFiltered_noPass_deepVar";
   my $gencode_bed_file = "/Shared/lss_tbraun/References/exons_gencode_33_20_merged.bed";
   my $dir;
   my $outfile = $file;
   $outfile =~ s/.*\/(.*)$/$1/;
   my $bgzipped_file = $outfile;
   $bgzipped_file = $outdir."/".$bgzipped_file.".gz";
   my $bed_filtered_file = $outdir."/bed_filtered_".$outfile;
   print "CHECK BGZIP: ".$bgzipped_file."\n";
   if (-e $bed_filtered_file) {
   } elsif (-e $bgzipped_file) {
      print "BGZIPPED FILE ALREADY EXISTS: ".$bgzipped_file."\n";
      my $check_csi = $bgzipped_file.".csi";
      print "CHECKING CSI FILE $check_csi\n"; 
      unless(-e $check_csi) {
         print "CSI FILE $check_csi DID NOT EXIST, SO CREATING1\n"; 
         system("bcftools index $bgzipped_file");
      }
   } else {
      my $to_bgzip = $bgzipped_file;
      $to_bgzip =~ s/\.gz$//;
      my $copy_cmd = "cp ".$file." ".$to_bgzip;
      print "80 ABOUT TO COPY: ".$copy_cmd."\n";
      system($copy_cmd);
      system("bgzip $to_bgzip");
      print "INDEX COMMAND: bcftools index $bgzipped_file\n";
      system("bcftools index $bgzipped_file");
   }
   my $name;
   my $cmd  ="";
    
   #$cmd = "bcftools view -R cut_exons_gencode33_20.bed $strelka_dir/$file -o strelka_filtered_by_bed_bcftools_gencode33_20/$outfile";
   print "BED_FILTERED_FILE $bed_filtered_file\n";
   if (-e $bed_filtered_file) {
      print "OUTPUT FILE ALREADY EXISTS\n";
   } else {
      $cmd = "bcftools view -R $gencode_bed_file $bgzipped_file -o $bed_filtered_file";
      print $cmd."\n";
      system($cmd);
   }
   my $IN;
   open $IN, $bed_filtered_file or die "Cannot open $bed_filtered_file :$!\n";
   my $line;
   my $past_headers = 0;
   my $line_count = 0;
   while ($line =<$IN>) {
      chomp $line;
      if ($past_headers > 0) {
         $line_count++;
         if ($line =~ m/PASS/){ 
            $count_pass++;
         }}
      if ($line =~ m/^#CHROM/) {
         $past_headers = 1
      } 

   }
   print "LINE COUNT: ".$line_count."\tCOUNT PASS: ".$count_pass."\n";
   return $line_count."\t".$count_pass;;
}
