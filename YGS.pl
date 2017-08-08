#!/usr/local/bin/perl

#   #  YGS.pl  based on k_mer_Solexa_bitvector_11b.pl (v_11b 8 Oct 2012 10AM)  A. Bernardo Carvalho  apr2007 v. 28dec2012 
use warnings;
use Bit::Vector;    #for use with kmer_size > 16 need to modify the source code (Toolbox.h)			
use Carp::Clan ;
use PerlIO::gzip;

											
												$program_version = "v_11b "."8 Oct 2012 10AM";
                                                                                        
#usage: perl YGS.pl kmer_size=15 mode=trace      trace=vir_1_f20.fasta.gz
#usage: perl YGS.pl kmer_size=15 mode=contig     contig=virCA.fasta.gz
#usage: perl YGS.pl kmer_size=15 mode=final_run  contig=virCA.fasta.gz  trace=fvir1f20m15k15rep5.vector_rep5.gz male_trace=virSangerq20m15k15.vector.gz  gen_rep=virCAm15k15.vector_rep2.gz


process_command_line("dummy");
#test_memory_use("dummy");
if ($kmer_size eq "not_set"){$kmer_size = 15 };   #default value: 15-mers.  user-set value: usually 18. Values above 16 require a modified Bit::Vector module, and possible a perl compiled with -Duse64bitint    perlbrew install  perl-5.14.2 -Duse64bitint
$bits =  4**$kmer_size  ;
if ($kmer_size > 16){if ($kmer_size > 19){die "\nexcessive kmer_size: $kmer_size. Normal values between 15 and 19\n"}
                    check_64bits_compatibility($kmer_size);
                    };
                    

					
if ($mode eq "contig" ){		$result_filename = $contig_prefix.".contig_result"; 
					$gen_rep_full_name = $contig_prefix.".gen_rep" ;
					$ctg_rep_full_name = $contig_prefix.".ctg_rep" ;
					($result_filename , $gen_rep_full_name , $ctg_rep_full_name)=avoid_overwriting($result_filename , $gen_rep_full_name , $ctg_rep_full_name);
					check_existence($contig_full_name);
					overture("dummy");
					populate_contig_vector("dummy");
					store_vector($vector_genome_repetitive,$gen_rep_full_name);
					store_vector($vector_btw_ctg_repetitive,$ctg_rep_full_name);
					};					
					
					
if ($mode eq "trace" ){$result_filename = $trace_prefix.".trace_result"; 
			check_existence($trace_fasta_name);
			$temp_name = $trace_prefix; $temp_name =~ s/_f\d+$//     ; $trace_raw_vector_name  = $temp_name .".trace"; #New in version 7 s_56_bCU3_f16 became s_56_bCU3
			$trace_filt_vector_name = $trace_prefix.".trace"; #New in version 7b:  will stay as s_56_bCU3_f16 
			$trace_raw_rep_vector_name  = $trace_raw_vector_name."_rep" ;      #New in version 7d
			$trace_filt_rep_vector_name = $trace_filt_vector_name."_rep" ;      #New in version 7d
			($result_filename , $trace_raw_vector_name , $trace_filt_vector_name , $trace_raw_rep_vector_name , $trace_filt_rep_vector_name )=avoid_overwriting($result_filename , $trace_raw_vector_name , $trace_filt_vector_name , $trace_raw_rep_vector_name , $trace_filt_rep_vector_name);
			overture("dummy");
# 			$softmasking = "false"; 					 #raw traces (ignores softmasking)
# 			populate_trace_vector_anysize("dummy"); 
# 			store_vector($vector_trace,$trace_raw_vector_name); 
# 			store_vector($vector_trace_rep,$trace_raw_rep_vector_name);
			$vector_trace  = Bit::Vector->new($bits);
			$vector_trace_rep = Bit::Vector->new($bits);
  			$softmasking = "true";  					 #filtered traces (implements softmasking)
			populate_trace_vector_anysize("dummy"); 
			store_vector($vector_trace,$trace_filt_vector_name);
#			store_vector($vector_trace_rep,$trace_filt_rep_vector_name);
			}; 
			


				
								


if ($mode eq "final_run" ){               
                if ($male_trace_full_name eq "not_used") {
                            $result_filename = $contig_prefix."_".$trace_prefix.".final_result";
                            ($result_filename)=avoid_overwriting($result_filename); 
                            check_existence($contig_full_name , $gen_rep_full_name ,  $trace_full_name);
                            #free_mem("before overture");
                            overture("dummy");
                            #free_mem("after  overture");
                            load_vector($genome_repetitive_loaded,$gen_rep_full_name);                            
                            #free_mem("after load_vector 1");
                            load_vector($vector_trace,$trace_full_name);   
                            #free_mem("after load_vector 2");
                            }
                else {
                            $result_filename = $contig_prefix."_".$trace_prefix."_".$male_trace_prefix.".final_result";
                            ($result_filename)=avoid_overwriting($result_filename);   
                            check_existence($contig_full_name , $gen_rep_full_name ,  $trace_full_name , $male_trace_full_name);
                             #free_mem("before overture");
                            overture("dummy");
                            #free_mem("after  overture");
                            load_vector($genome_repetitive_loaded,$gen_rep_full_name);
                            load_vector($vector_trace,$trace_full_name);
                            load_vector($vector_male_trace,$male_trace_full_name);
                            }
                populate_contig_vector("dummy"); 
                if ($save_memory ne "yes"){genome_wide_analysis("dummy")};      #  for running k18 use option save_memory=yes . Saves memory by  not using the cummulative bit-array used in genome_wide_analysis
                };
                

				
print "\n\nprogram finished at: ",scalar localtime,"\n\n\n";
close (RESULT);
exit;

#read(STDIN,$i,1);




sub check_existence{
my @file = @_ ;
$flag_input_file_error = "false";
#print "check_existence @file";
for ($i = 0 ; $i < @file ; $i++){
#print "file $i   $file[$i]\n\n";
	if  (not( -e $file[$i] )) {
		print "\nfile $file[$i] do not  exist. Please check name.\n";
		$flag_input_file_error = $file[$i];
		};
	};
if ($flag_input_file_error ne "false"){print "\n\n input file error. File $flag_input_file_error does not exist.\n\n"; die};
return};



sub avoid_overwriting{					# avoids overwriting files ?check if input files exist? 
my @file = @_ ;
#print "@file\n";
#my $filename="";
for ($i = 0 ; $i < @file ; $i++){
	$flag_overwrite = "false";
#print "current file: $file[$i]\n";
	while ( ( -e $file[$i] ) && ($flag_overwrite eq "false")) {
#print "current file: $file[$i]\n";
		print "\nfile $file[$i] already exists. Overwrite ? (y/n)\n";
		read(STDIN,$answer,1); chomp($answer); while ($answer eq ""){read(STDIN,$answer,1)};
		if ($answer eq "n"){
			print "\ninput new filename:\n"; 
			$file[$i] = readline(STDIN); chomp($file[$i]);
			while ($file[$i] eq ""){$file[$i] = readline(STDIN); chomp($file[$i]);} 
			print "filename: $file[$i]\n" ;
			}
		if ($answer eq "y"){
			print "\nfile $file[$i] will be overwrited\n";
			$flag_overwrite="true";
			}
		}
	}
return @file}


sub process_command_line{
$mode = ''; $contig_prefix = ''; $contig_full_name = ''; $trace_prefix = ''; $trace_full_name = ''; $gen_rep_full_name = ''; $ctg_rep_full_name = '';
$trace_raw_vector_name = ''; $trace_filt_vector_name = ''; $trace_fasta_name = ''; $kmer_size="not_set"; $male_trace_full_name="not_used"; $male_trace_prefix=""; $save_memory="no";
if (@ARGV == 0){usage_v_8("dummy"); die}
$num_arg = @ARGV;
for ($i = 0 ; $i < $num_arg ; $i++){
	#print "\n\n num_arg: $num_arg     command argument: @ARGV    \nargument number $i  $ARGV[$i]\n\n";	
	$arg = $ARGV[$i];
	my($flag_match)=0;
	if ($arg =~ m/mode=(\w+)/ ) { $mode = $1; $flag_match ++};  
	if ($arg =~ m/contig=(\w+)/ ){ $contig_prefix = $1; $flag_match ++};  
	if ($arg =~ m/contig=([\w\.]+)/ ){$contig_full_name = $1; $flag_match ++}; 	 
	if ($arg =~ m/\btrace=(\w+)/ ){ $trace_prefix = $1; $flag_match ++};    #\b is word boundary. This is to avoid matching the v.11 male_trace=
	if ($arg =~ m/\btrace=([\w\.]+)/ ){$trace_fasta_name = $1; $flag_match ++;
					$trace_full_name = $trace_fasta_name};
	if ($arg =~ m/gen_rep=([\w\.]+)/ ){ $gen_rep_full_name = $1; $flag_match ++ };      
    if ($arg =~ m/kmer_size=(\d+)/ ){ $kmer_size = $1; $flag_match ++ };
    if ($arg =~ m/\bmale_trace=(\w+)/ ){$male_trace_prefix = $1; $flag_match ++};    
    if ($arg =~ m/\bmale_trace=([\w\.]+)/ ){$male_trace_full_name = $1; $flag_match ++};
    if ($arg =~ m/save_memory=(\w+)/ ) { $save_memory = $1; $flag_match ++;
                    if ( ($save_memory ne "yes")&&($save_memory ne "no")  ) {die "\nUse either save_memory=yes or save_memory=no. \n"}                                        
                    };
    if ($flag_match==0){die "Unrecognized option $arg . Check spelling. \n"};
	};
#print "mode: $mode   contig_prefix: $contig_prefix  contig_full_name: $contig_full_name trace_prefix: $trace_prefix trace_full_name: $trace_full_name gen_rep_full_name: $gen_rep_full_name ctg_rep_full_name: $ctg_rep_full_name hash_full_name: $hash_full_name cut_off: $cut_off\n\n";
if ($mode eq "contig"){
	if ( ($contig_prefix eq "") || ($trace_prefix ne "") || ($gen_rep_full_name ne "")  ){usage_v_8("dummy"); die "mode $mode is incompatible with choosen files. Check usage above\n"}
	}
if ($mode eq "trace"){
	if ( ($contig_prefix ne "") || ($trace_prefix eq "") || ($gen_rep_full_name ne "")  ){usage_v_8("dummy"); die "mode $mode is incompatible with choosen files. Check usage above\n"}
	}
if ($mode eq "final_run"){
	if ( ($contig_prefix eq "") || ($trace_prefix eq "") || ($gen_rep_full_name eq "")  ){usage_v_8("dummy"); die "mode $mode is incompatible with choosen files. Check usage above\n"}
	}

return}


sub overture{ 
open (RESULT, ">$result_filename") || die "not possible to create the file $result_filename\n";
if ($mode ne "interactive" ){select RESULT; $| = 1};  
print "\n\n",scalar localtime,"\n\n";
if ($mode eq "contig"){print  "Input\ncontig file (fasta): $contig_full_name\n\nOutput\ngenome repetitive file (vector): $gen_rep_full_name  \ncontig repetitive file (vector): $ctg_rep_full_name    \npartial results (text): $result_filename\n\n"};
if ($mode eq "trace"){print  "Input\ntrace file (fasta): $trace_fasta_name \n\nOutput\ntrace files (vector):   $trace_raw_vector_name          $trace_filt_vector_name\npartial results (text): $result_filename\n\n"};
if ($mode eq "final_run"){print  "Input\ncontig file (fasta): $contig_full_name     \ntrace file (vector): $trace_full_name  \ngenome repetitive file (vector): $gen_rep_full_name  \nmale trace file (vector): $male_trace_full_name \n\nOutput \nfinal results (text): $result_filename\n\n"};


#print "contig file: $contig_full_name \nhash_repetitive file:  $gen_rep_full_name \ntrace_file: $trace_file \n\n";
print "program: $0  version: $program_version    PID:  $$    program mode: $mode \nkmer_size: $kmer_size\n\n";
print "command line:\n\nperl $0 @ARGV\n\n";
if ( ($mode eq "hash_repetitive")||($mode eq "contig" ) ){
$vector_genome_repetitive = Bit::Vector->new($bits);
$vector_btw_ctg_repetitive = Bit::Vector->new($bits);  #needed for saving ctg_rep
$vector_genome = Bit::Vector->new($bits); 
}
if ($mode eq "trace" ){
	$vector_trace  = Bit::Vector->new($bits);
	$vector_trace_rep = Bit::Vector->new($bits);   
	};
if ($mode eq "final_run"){
	$genome_repetitive_loaded = Bit::Vector->new($bits);
	$vector_trace = Bit::Vector->new($bits);
	if ($male_trace_full_name ne "not_used"){$vector_male_trace= Bit::Vector->new($bits) };  # if male traces are not used, this saves the memory
    if ($save_memory ne "yes"){$vector_genome = Bit::Vector->new($bits)};   #this bit-array is used only for the genome_wise_analysis procedure 
	};
return};



sub populate_contig_vector {						#Read contig file and populate contig vector with k-mers	
$contig_num = 0;
$contig_fasta = '';	
$fasta_size = 0; 
$total_kmer = 0;
$total_max_kmer = 0;
if ($mode ne "interactive" ){open DTF, "<:gzip(autopop)", "$contig_full_name"  || die "not possible to open the file $contig_full_name\n"}
			else{open DTF, "<:gzip(autopop)", "$curr_ctg" || die "not possible to open the file $curr_ctg\n"};
while (<DTF>) {
	chomp($_);
	if ($_ =~ /^>/){
		if ($contig_num == 0){						#first contig of the file
			$defline = $_ ;
			$contig_num ++ ;
			}
        else{					                    #new contig found. process the old, and update the contig name, etc
            process_contig_fasta ($contig_fasta);                    
            $contig_num ++;
            $contig_fasta = '';
            };
		};	
	if ($_ !~ /^>/){
		$curr_line = $_ ; 
		if (not( ($mode eq "interactive" )||($mode eq "mask_fasta" ) || ($mode eq "mask_fasta2" ) ) ){$curr_line =~ s/N+/n/gi};		#substitute lines of pure NNNNNNNNNNNNNNNN by one N
		$contig_fasta = $contig_fasta.$curr_line ;
		};
}						

process_contig_fasta ($contig_fasta) ;                  #process the last contig    
print "\n\nprocessing contig finished \ncontig_wide analysis: (Note: kmers_found is the sum of kmers found in each contig, and hence contains between-contig repeats)\n\n$contig_num  contigs.      $total_max_kmer  max_kmers      $total_kmer kmers_found";
if (($mode eq "hash_repetitive" )||($mode eq "contig" )){
	$kmer_repetitive = $vector_genome_repetitive->Norm();
	print "        $kmer_repetitive  repetitive_kmers_found\n\n";
	}else {print "\n\n"};
print "within contigs repeats: ", ($total_max_kmer - $total_kmer) , " kmers\n\n";
return};
	


sub process_contig_fasta {
no warnings 'portable';   #avoids a huge number of warnings in the error file (at UW cluster): Binary number > 0b11111111111111111111111111111111 non-portable 
$defline =~ m/(^\S+)/; $defline_short = $1;    # m// implicitely use the  $_ variable_  
my ($fasta) = @_;
$ctg_max_kmer = 0;
$fasta_size = length($fasta); 
$vector_contig = Bit::Vector->new($bits);   #new in v10
#print "$defline \n    $fasta_size \n"; 


for ($i = 0 ; ($i < $fasta_size - ($kmer_size-1) ) ; $i++){
	$curr_15 = substr($fasta, $i , $kmer_size) ;            #I mantained the old variable name, but from v 9 onwards  kmer size can be different from 15
	#print "$curr_15\n";	
	chomp($curr_15);
	$seq_or = uc($curr_15) ; 
		if (not($seq_or =~ m/[^ATGC]/)){				     #skips kmer containing anything different from ATGC (e.g., N, wrong characters, etc)
			$ctg_max_kmer ++ ;
			$seq_rc = reverse $seq_or ;
			$seq_rc =~ tr/ATCG/TAGC/;
			if ($seq_or lt $seq_rc){$s15bin = $seq_or}
			else{$s15bin = $seq_rc};
			$s15bin =~ s/A/00/g ; $s15bin =~ s/T/11/g ; $s15bin =~ s/G/10/g ; $s15bin =~ s/C/01/g ;
			$s15dec = oct("0b" . $s15bin);

            if ($mode ne "final_run" ){                             #final_run is the slowest mode. this avoids checking all mode options  , for each kmer


                    if ($mode eq "contig"){
                        if ( $vector_genome->bit_test($s15dec) ){
                            $vector_btw_ctg_repetitive->Bit_On($s15dec);
                            $vector_genome_repetitive->Bit_On($s15dec);
                            }
                        if ( $vector_contig->bit_test($s15dec) ) {
                            $vector_genome_repetitive->Bit_On($s15dec)
                            }				
                        }


                    }
			
			$vector_contig->Bit_On($s15dec);
			}else{}
			#}else{print "non ATGC bases detected in contig $defline : $seq_or\n" }
	}

	
	

if ( ($save_memory eq "no")||($mode ne "final_run") ) { $vector_genome->Or($vector_genome,$vector_contig)	}			#Add the current contig data to the genome vector (cumulative contig)
$ctg_kmer = $vector_contig->Norm();						#number of set bits in the vector. In my case, the number of different k-mer found; 
$total_kmer = 	$total_kmer + $ctg_kmer;
$total_max_kmer = $total_max_kmer + $ctg_max_kmer;

if ($mode eq "final_run" ){             #modified in v11a (to allow use of male traces) and v11b (to reduce memory use)
    $vector_ctg_unmatched = Bit::Vector->new($bits);
    $vector_ctg_unmatched->AndNot($vector_contig,$vector_trace); 
    $ctg_unmatched_kmer = $vector_ctg_unmatched->Norm();
    undef($vector_ctg_unmatched);
    $vector_ctg_scp = Bit::Vector->new($bits);
    $vector_ctg_scp->AndNot($vector_contig,$genome_repetitive_loaded) ; $ctg_sc_kmer = $vector_ctg_scp->Norm();
    undef($vector_contig);  #new in v10, to save memory
    if ($male_trace_full_name eq "not_used"){           #validating male traces NOT used      
            $vector_ctg_scp_unmatched = Bit::Vector->new($bits);
            $vector_ctg_scp_unmatched->AndNot($vector_ctg_scp,$vector_trace); $ctg_sc_unmatched_kmer = $vector_ctg_scp_unmatched->Norm();    
            undef($vector_ctg_scp);
            undef($vector_ctg_scp_unmatched);
            if ($contig_num == 1){print "GI                  NUM    MAX_K    K    UK    SC_K    SC_UK    P_SC_UK \n"};           #names simplified in v11
            if ($ctg_sc_kmer > 0){ $P_SC_UK = sprintf("%.1f" , 100*$ctg_sc_unmatched_kmer / $ctg_sc_kmer) }else{$P_SC_UK="."};   #new in v10
            print "$defline_short $contig_num  $ctg_max_kmer $ctg_kmer $ctg_unmatched_kmer $ctg_sc_kmer  $ctg_sc_unmatched_kmer  $P_SC_UK \n";
            }
    else{                                       #validating male traces used         
            #free_mem("before vector_validated_ctg_scp");
            $vector_validated_ctg_scp = Bit::Vector->new($bits);
            $vector_validated_ctg_scp->And($vector_ctg_scp,$vector_male_trace);
            #free_mem("after  vector_validated_ctg_scp");
            $validated_ctg_sc_kmer = $vector_validated_ctg_scp->Norm(); 
            undef($vector_validated_ctg_scp);
            $vector_ctg_scp_unmatched = Bit::Vector->new($bits);
            $vector_ctg_scp_unmatched->AndNot($vector_ctg_scp,$vector_trace); 
            undef($vector_ctg_scp);
            $ctg_sc_unmatched_kmer = $vector_ctg_scp_unmatched->Norm();
            $vector_validated_ctg_scp_unmatched = Bit::Vector->new($bits);
            #free_mem("after  vector_validated_ctg_scp_unmatched");
            $vector_validated_ctg_scp_unmatched->And($vector_ctg_scp_unmatched,$vector_male_trace);
            undef($vector_ctg_scp_unmatched);
            $validated_ctg_sc_unmatched_kmer = $vector_validated_ctg_scp_unmatched->Norm();
            undef($vector_validated_ctg_scp_unmatched);
            if ($contig_num == 1){print "GI                  NUM    MAX_K    K    UK    SC_K    SC_UK    P_SC_UK   VSC_K   VSC_UK   P_VSC_UK \n"}; 
            if ($ctg_sc_kmer > 0){ $P_SC_UK = sprintf("%.1f" , 100*$ctg_sc_unmatched_kmer / $ctg_sc_kmer) }else{$P_SC_UK="."};
            if ($validated_ctg_sc_kmer > 0){ $P_VSC_UK = sprintf("%.1f" , 100*$validated_ctg_sc_unmatched_kmer / $validated_ctg_sc_kmer) }else{$P_VSC_UK="."};            
            print "$defline_short $contig_num  $ctg_max_kmer $ctg_kmer $ctg_unmatched_kmer $ctg_sc_kmer  $ctg_sc_unmatched_kmer  $P_SC_UK  $validated_ctg_sc_kmer   $validated_ctg_sc_unmatched_kmer  $P_VSC_UK \n";
            }
    }

if (($mode eq "hash_repetitive" )||($mode eq "contig" )){
	if ($contig_num == 1){print "defline_short contig_num  contig_max_kmer contig_kmer\n"};
	print "$defline_short $contig_num  $ctg_max_kmer $ctg_kmer\n";
	};

	

    
$defline = $_ ;
#$contig_fasta = ''; $fasta = "";
undef($contig_fasta) ; undef($fasta);  #new in v 11b
return;
}



sub populate_trace_vector_anysize{			#Read trace file  and populate trace_vector with k-mers.	
my $trace_num = 0;	
my $trace_size = 0; 
my $trace_size_total = 0;
my $trace_valid_kmers_total = 0;  
#open (DTF, "< $trace_fasta_name") || die "not possible to open the file $trace_fasta_name\n";
open DTF, "<:gzip(autopop)","$trace_fasta_name" or die "not possible to open the file $trace_fasta_name\n";
while (<DTF>) {
	if ($_ =~ /^>/){$_ = readline(DTF);}					#skips the deflines
	$trace_num ++ ;
#print STDOUT $_ ;	
	chomp($_);
	if ($softmasking eq "false"){	$_ = uc($_)};				#NEW in version 7: inhibits soft masking.
	$trace_size = length($_); $trace_size_total = $trace_size_total + $trace_size;
	for ($i = 0 ; ($i < $trace_size -($kmer_size-1) ) ; $i++){
		$curr_15 = substr($_, $i , $kmer_size) ;
		#print "$curr_15\n";	
		chomp($curr_15);					
		$seq_or =    ($curr_15) ; 
#print "$softmasking $_  $seq_or  ";					
		if (not($seq_or =~ m/[^ATGC]/)){				#skips 15-mer containing anything different from ATGC (e.g., N, wrong characters, etc)
			$trace_valid_kmers_total ++ ;  	#new in v. 7d: counts the total number of valid kmers (many are identical)
			$seq_rc = reverse $seq_or ;
			$seq_rc =~ tr/ATCG/TAGC/;
			if ($seq_or lt $seq_rc){$s15bin = $seq_or}
			else{$s15bin = $seq_rc};
			$s15bin =~ s/A/00/g ; $s15bin =~ s/T/11/g ; $s15bin =~ s/G/10/g ; $s15bin =~ s/C/01/g ;
			$s15dec = oct("0b" . $s15bin); 
#print "$s15dec   $seq_or   \n ";
			if ( $vector_trace->bit_test($s15dec) ) {		#implements repetitive traces (those appearing more than once in the trace file. To avid spurious hits due to sequencing errors 
				$vector_trace_rep->Bit_On($s15dec);
				}				
			$vector_trace->Bit_On($s15dec);
			}else{}
			#}else{print "non ATGC bases detected in trace $trace_num : $seq_or\n" }
		}
	};
$trace_kmers_number = $vector_trace->Norm();
print "\n\nprocessing traces finished   file: $trace_fasta_name .  softmasking mode: $softmasking  $trace_num traces.    $trace_size_total total raw bp    $trace_valid_kmers_total total valid kmers (many will be identical)     $trace_kmers_number trace_kmers (not trace sequences!) found\n";
close (DTF);
return};







sub load_vector{			#two parameters: vector to be loaded, and filename
($vector,$vector_filename) = @_;
my $kmer_number;
open VECTOR_L, "<:gzip(autopop)","$vector_filename" or die "not possible to open the file $vector_filename\n";
$string = readline(VECTOR_L);
$vector->from_Hex($string);
$kmer_number = $vector->Norm();
#$string = "";								#to save memory!
undef($string);         #new in v11b  . This really saves memory
print "\n\nload_vector finished. $kmer_number kmers loaded from file $vector_filename  \n\n";
if ($kmer_number == 0){print "\n Zero kmers found in file $vector_filename . This may be caused by your perl installation not supporting 64-bits. Try  perlbrew install  perl-5.16.0 -Duse64bitint . Or reduce kmer size. \n"};
my($vector_size)= $vector->Size();
my($predicted_vector_size) = 4**$kmer_size ;  #changed in 9c. 
if ($vector_size != $predicted_vector_size ){
                            print "\nvector size in file $vector_filename ( $vector_size ) differs from the vector size  expected ($predicted_vector_size) from the declared kmer_size( $kmer_size ) \n";
                            my($inferred_kmer_size) = log($vector_size) / log(4) ; 
                            print "vector size in file $vector_filename ( $vector_size ) is compatible with a kmer size of $inferred_kmer_size\n"; 
                            die;
                            };
close (VECTOR_L);
#undef($vector); undef(@_);  # failed to reduce memory usage (as judged by free)
return};






sub store_vector{        		#two parameters: vector to be stored, and filename
($vector,$vector_filename) = @_;
my $kmer_number;
#open (VECTOR_S, ">$vector_filename") || die "not possible to create the file $vector_filename\n";
open VECTOR_S, ">:gzip","$vector_filename\.gz" or die "not possible to create the file $vector_filename\n";
$string = $vector->to_Hex();
print VECTOR_S $string;
#$string = "";								#to save memory!
undef($string);     #new in v 11b 8oct
$kmer_number = $vector->Norm();
print "\n\nstore_vector finished. $kmer_number kmers stored in  file $vector_filename\n";
close (VECTOR_S);
$vector->Resize(1);
return};







sub genome_wide_analysis{				
#$vector_ctg_scp->Resize(1); $vector_ctg_unmatched->Resize(1);$vector_ctg_scp_unmatched->Resize(1); #to save memory , I am resizing (to 1 bit) the no longer needed vectors
undef($vector_ctg_scp); undef($vector_ctg_unmatched); undef($vector_ctg_scp_unmatched);
$vector_temp = Bit::Vector->new($bits);
$vector_genome_scp = Bit::Vector->new($bits);

$genome_kmer = $vector_genome->Norm();     										#total kmers in the genome
$vector_genome_scp->AndNot($vector_genome,$genome_repetitive_loaded) ; $genome_sc_kmer = $vector_genome_scp->Norm();     	#sc kmers in the genome
$genome_rep_kmer = $genome_repetitive_loaded->Norm();									#repetitive kmers in the genome
$vector_temp->AndNot($vector_genome,$vector_trace); $genome_unmatched_kmer = $vector_temp->Norm();	# unmatched total kmers in the genome

$vector_temp->AndNot($vector_genome_scp,$vector_trace); $genome_unmatched_sc_kmer = $vector_temp->Norm();#unmatched sc  kmers in the genome
$vector_temp->AndNot($genome_repetitive_loaded,$vector_trace); $genome_unmatched_rep_kmer = $vector_temp->Norm();#unmatched rep kmers in the genome

$vector_temp->AndNot($vector_trace,$vector_genome); $trace_unmatched_kmer = $vector_temp->Norm();	# unmatched kmers in the traces
$trace_total_kmer = $vector_trace->Norm();
print "\n\ngenome_wide analysis:\n\n";
print "total kmers in the genome: $genome_kmer\n";
print "scp kmers in the genome:   $genome_sc_kmer\n";
print "rep kmers in the genome:   $genome_rep_kmer"; printf("%s%.1f%s" , " ( ", 100*$genome_rep_kmer/$genome_kmer , " % of the total )\n\n"); 
print "unmatched total kmers in the genome: $genome_unmatched_kmer"; printf("%s%.1f%s" , " ( ", 100*$genome_unmatched_kmer/$genome_kmer , " % of the total genome kmers)\n");
print "unmatched scp kmers in the genome:   $genome_unmatched_sc_kmer"; printf("%s%.1f%s" , " ( ", 100*$genome_unmatched_sc_kmer/$genome_sc_kmer , " % of the scp genome kmers )\n");
print "unmatched rep kmers in the genome:   $genome_unmatched_rep_kmer"; printf("%s%.1f%s" , " ( ", 100*$genome_unmatched_rep_kmer/$genome_rep_kmer , " % of the rep genome kmers )\n\n");
print "total kmers in the traces:     $trace_total_kmer\n";
print "unmatched kmers in the traces: $trace_unmatched_kmer";  printf("%s%.1f%s" , " ( ", 100*$trace_unmatched_kmer/$trace_total_kmer , " % of the total )\n"); 
undef($vector_temp);
return};



sub usage_v_11{
print "\n\nusage: perl YGS.pl mode=contig/trace/final_run contig=contig_file   gen_rep=genome_repetitive_file  trace=trace_file   \n\nExamples:\n";
print "perl YGS.pl mode=contig  contig=wgs3_all.txt \n\n";
print "perl YGS.pl mode=trace  trace=s_1.fasta \n\n";
print "perl YGS.pl mode=final_run  contig=wgs3_all.txt   trace=s_1.trace   gen_rep=wgs3_all.gen_rep  \n\n";
print "\n\n\n\nThree modes of run:\n\n\n";
print "contig mode:           Input: contig fasta file\n                       Output: *.gen_rep.gz file, *.ctg_rep.gz file  (and *.contig_result file)\n\n";
print "trace mode:            Input: trace fasta  file\n                       Output: *.trace.gz file  ,  *.trace_rep.gz file  (and .trace_result file)\n\n";
print "final_run  mode:       Input: contig fasta file,  *.trace file, *.gen_rep file, *.ctg_rep file\n                   Output: *.final_result file\n\n\n";
return};






sub check_64bits_compatibility{   #kmer_size bigger than 16 requires modified Bit::Vector and probably perl compilation with -Duse64bitint  (e.g. ,  perlbrew install  perl-5.14.2 -Duse64bitint)
no warnings 'portable';   #avoids one error message: Binary number > 0b11111111111111111111111111111111 non-portable 
my ($k_size) = @_;
my($flag_64bits_incompatibility)=0;
my($check_integer) = 2**36 ;                                    #checks perl exponentiation above 2^32
if ($check_integer ne "68719476736"){
        print "ERROR: this perl installation  fails in exponentiation above 2^32\n";
        $flag_64bits_incompatibility ++ ;
        }
        
print "passed  exponentiation above 2^32\n";
 
my($check_binary) = "10000000000000000000000000000000000" ;     # checks perl bin->dec conversion above 2^32). The binary number is 2^34
$check_integer = oct("0b" . $check_binary);
if ($check_integer ne "17179869184"){
        print "ERROR: this perl installation fails in bin->dec conversion above 2^32 \n";
        $flag_64bits_incompatibility ++ ;
        }
else{print "passed bin->dec conversion above 2^32 \n"};        
        
my($Long_Bits) = Bit::Vector->Long_Bits() ;                     #checks installed Bit::Vector module. 
my($Word_Bits) = Bit::Vector->Word_Bits() ;
if ( ($Long_Bits != 64) || ($Long_Bits != 64)  ){
        print "ERROR: installed Bit::Vector module does not support kmer_size above 16 (Long_Bits / Word_Bits test)\n";
        $flag_64bits_incompatibility ++ ;
        die "Not possible to run the program in this system due to large kmer_size ($k_size).\nkmer_size bigger than 15 requires modified Bit::Vector and probably perl compilation with -Duse64bitint  (e.g. ,  perlbrew install  perl-5.14.2 -Duse64bitint)\n";
        }
else{print "passed    Long_Bits / Word_Bits=64\n"};      

my($bits) =  4**$k_size  ;
$check_vector = Bit::Vector->new($bits);                   #functional tests of installed Bit::Vector module.
my($vector_size) = $check_vector->Size();
if ($bits != $vector_size){
        print "ERROR: installed Bit::Vector module does not support kmer_size above 15. failed in new() or Size()\n";
        $flag_64bits_incompatibility ++ ;
        }
else{print "passed    functional test of installed Bit::Vector module  new()  Size()\n"};         
        
my($state_before_On,$state_after_On ) = 2;                      #functional tests of installed Bit::Vector module.
$state_before_On = $check_vector->bit_test($vector_size -1);
$check_vector->Bit_On($vector_size -1);
$state_after_On  = $check_vector->bit_test($vector_size -1);
if ( ($state_before_On != 0) || ($state_after_On != 1)  ){
        print "ERROR: installed Bit::Vector module does not support kmer_size above 16. failed in Bit_On() or  bit_test() \n";
        $flag_64bits_incompatibility ++ ;
        }   
else{print "passed    functional test of installed Bit::Vector module  bit_test()  Bit_On()\n"};
 
if ($flag_64bits_incompatibility > 0){die "Not possible to run the program in this system due to large kmer_size ($k_size).\nkmer_size bigger than 16 requires modified Bit::Vector and probably perl compilation with -Duse64bitint  (e.g. ,  perlbrew install  perl-5.16.0 -Duse64bitint)\n"}
else{print "passed all 64-bits tests\n"};
undef($check_vector);  #to save memory
return}




