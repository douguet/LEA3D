#!/usr/bin/perl


	local($fileopen)=@ARGV;

	$leaexe=$0;
	#windows requires 2 steps
        $leaexe=~s/lea3d-MAKE_FGTS\.pl$//;
	$leaexe=~s/\/$//;
        #print "perl scripts in $leaexe\n";

	if($fileopen eq ''){
		die "file ?\nusage: makefgt <file.sdf> (make_fgts.sdf = ring + fused_rings + special + linker + substituent ; acyclic.sdf=substituent.sdf+linker.sdf; if 1 fgt: 1 acyclic then fills acyclic.sdf; if 1 ring then fills ring file and make_fgts.sdf)\n";
	};

$cyclesdf=$leaexe."/lea3d-CYCLE_SDF.pl";
require $cyclesdf;

$flagnew=1;
$flagmend=0;
$nomol=0;

unlink "ring.sdf" if(-e "ring.sdf");
unlink "fused_rings.sdf" if(-e "fused_rings.sdf");
unlink "special.sdf" if(-e "special.sdf");
unlink "substituent.sdf" if(-e "substituent.sdf");
unlink "acyclic.sdf" if(-e "acyclic.sdf");
unlink "linker.sdf" if(-e "linker.sdf");

unlink "linker_tmp.sdf";
unlink "cycle_tmp.sdf";
unlink "cyclefus_tmp.sdf";
unlink "special_tmp.sdf";
									

$fcycle="ring.sdf";
$fcyclefus="fused_rings.sdf";
$flinker="acyclic.sdf";
$fspecial="special.sdf";

#$fcyclelinker="rings_and_linker.sdf";
$fsubstituent="substituent.sdf";
$ftruelinker="linker.sdf";

open(IN,"<$fileopen");
while (<IN>){

	$conv2=$_;

	if ($flagnew){
		$compt=0;
		$i=1;
		$j=0;
		$blanc=" ";
		@bond='';
		@listb='';
		
		@ifonc='';
		
		@typeb='';
		@type='';
		@corx='';
		@cory='';
		@corz='';
		@lignecor='';
		@lignebond='';
		$lignedata='';
		$flagdata=0;
		$flagnew=0;
		$differentiel=0;
		$differentielx=0;
		$nomol++;
	};

	@get = split(' ',$_);
	$compt++;

	
	if (($compt > 4) && ($i <= $nbatom)){
		$corx[$i]=$get[0];
		$cory[$i]=$get[1];
		$corz[$i]=$get[2];
		$type[$i]=$get[3]; 
		
		if($type[$i] eq 'X'){
			print "Warning \'X\' dummy atom detected and changed to H !\n";
			$type[$i]="H";
			$_=~s/ X / H /;
		};	

		$lignecor[$i]=$_;
		$i++;
	};

	if (($compt > 4) && ($i > $nbatom) && ($j <= $nbond)){
		if ($j == 0){
			$j++;
		}
		else{
		        	@coller=split(' *',$get[0]);
				@coller2=split(' *',$get[1]);
				if(@coller==6 && $get[1] ne ""){
					$get[0]=$coller[0].$coller[1].$coller[2];
					$get[2]=$get[1];
					$get[1]=$coller[3].$coller[4].$coller[5];
				}
				elsif(@coller==6 && $get[1] eq ""){
					$get[0]=$coller[0].$coller[1];
					$get[1]=$coller[2].$coller[3].$coller[4];
					$get[2]=$coller[5];
				}
				elsif(@coller==5){
					if($_=~/^\s/){
						$get[0]=$coller[0].$coller[1];
						$get[2]=$get[1];
						$get[1]=$coller[2].$coller[3].$coller[4];
					}
					else{
						$get[0]=$coller[0].$coller[1].$coller[2];
						$get[2]=$get[1];
						$get[1]=$coller[3].$coller[4];
					};
                                }
                                elsif(@coller==4){
                                        if($_=~/^\s/){
                                                $get[0]=$coller[0];
                                                $get[2]=$get[1];
                                                $get[1]=$coller[1].$coller[2].$coller[3];
                                        }
                                        else{
                                                $get[0]=$coller[0].$coller[1].$coller[2];
						$get[2]=$get[1];
                                                $get[1]=$coller[3];
                                        };					
				}
				elsif(@coller2==4){
					$get[1]=$coller2[0].$coller2[1].$coller2[2];
					$get[2]=$coller2[3];
				}
				elsif(@coller==7){
					$get[0]=$coller[0].$coller[1].$coller[2];
					$get[1]=$coller[3].$coller[4].$coller[5];
					$get[2]=$coller[6];
				};

			$bond[$get[0]]=$bond[$get[0]].$blanc.$get[1].$blanc.$get[2];
			$listb[$get[0]]=$listb[$get[0]].$blanc.$get[1];
			$typeb[$get[0]]=$typeb[$get[0]].$blanc.$get[2];

			$bond[$get[1]]=$bond[$get[1]].$blanc.$get[0].$blanc.$get[2];
			$listb[$get[1]]=$listb[$get[1]].$blanc.$get[0];
			$typeb[$get[1]]=$typeb[$get[1]].$blanc.$get[2];

			$ifonc[$get[0]]=$ifonc[$get[0]].$blanc.$get[1].$blanc;
			$ifonc[$get[1]]=$ifonc[$get[1]].$blanc.$get[0].$blanc;
			
			$lignebond[$j]=$_;
			$j++;
		};
	};
	
	if($compt > 4 && $i > $nbatom && $j > $nbond){
		$flagdata=1 if($_=~/^>/ && ($_=~/CAS/ || $_=~/cas/ || $_=~/NAME/ || $_=~/MDLNUMBER/ || $_=~/ZINC/));
	};	
	
	if ($compt == 4){
		$nbatom=$get[0];
		$nbond=$get[1];

			@coller=split(' *',$nbatom);
			if(@coller > 3 && @coller == 6){
				$nbatom=$coller[0].$coller[1].$coller[2];
				$nbond=$coller[3].$coller[4].$coller[5];
			}
			elsif( @coller > 3 && @coller == 5){
				if($_=~/^\s/){
					$nbatom=$coller[0].$coller[1];
					$nbond=$coller[2].$coller[3].$coller[4];
				}
				else{
					$nbatom=$coller[0].$coller[1].$coller[2];
					$nbond=$coller[3].$coller[4];
				};
			};

		#$compt++;
	};

	if ($conv2 =~/\$\$\$\$/){
		$flagdata=0;
		$flagnew=1;
		&sdf;
		$flagmend=0;
		$differentield=$nbond - $differentiel+($differentielx /2);

		if($differentield == 0){
			open(OUTC,"<cycle_tmp.sdf");
			open(OUTD,">>$fcycle");
			while(<OUTC>){
				print OUTD $_;
			};
			close(OUTC);
			close(OUTD);
			unlink "cycle_tmp.sdf";

                        open(OUTC,"<cyclefus_tmp.sdf");
                        open(OUTD,">>$fcyclefus");
                        while(<OUTC>){
                                print OUTD $_;
                        };
                        close(OUTC);
                        close(OUTD);
			unlink "cyclefus_tmp.sdf";

			open(OUTC,"<special_tmp.sdf");
			open(OUTD,">>$fspecial");
			while(<OUTC>){
				print OUTD $_;
			};
			close(OUTC);
			close(OUTD);
			unlink "special_tmp.sdf";	
			
                        open(OUTC,"<linker_tmp.sdf");
                        open(OUTD,">>$flinker");
			$nbx=0;
			$nbxi=0;
			@nbxtab='';
                        while(<OUTC>){
                                if($_=~/\$\$\$\$/){
					$nbxtab[$nbxi]=$nbx;
					$nbx=0;
					$nbxi++;
				};
				print OUTD $_;
				@getx = split(' ',$_);
				$nbx++ if($getx[3] eq "X");
                        };
                        close(OUTC);
                        close(OUTD);
			$nbxi=0;
			open(OUTC,"<linker_tmp.sdf");
			open(OUTF,">>$fsubstituent");
			open(OUTG,">>$ftruelinker");
			while(<OUTC>){
				if($nbxtab[$nbxi]==1){
					print OUTF $_;
				}
				elsif($nbxtab[$nbxi] > 1){
					print OUTG $_;
				};
				if($_=~/\$\$\$\$/){
					$nbxi++;
				};
			};
			close(OUTC);
			close(OUTF);
			close(OUTG);
			unlink "linker_tmp.sdf";
		}
		else{
			#print "Molecule $nomol excluded because (cycle > 10 atoms ?) PROBLEM nbbond read = $nbond ==  $differentiel +  $differentielx /2 \n";
			unlink "linker_tmp.sdf";
			unlink "cycle_tmp.sdf";
			unlink "cyclefus_tmp.sdf";
			unlink "special_tmp.sdf";
		};
	}
	else{
		$flagmend=1;
	};
	if($flagdata){
		$lignedata=$lignedata.$_;
	};

};
close(IN);

#create make_fgts.sdf with all fragments
unlink "make_fgts.sdf" if(-e "make_fgts.sdf");
open(UYT,">make_fgts.sdf");
#system("touch make_fgts.sdf");
if(-e "ring.sdf" && !-z "ring.sdf"){
	open(DFH,"<ring.sdf");
	while(<DFH>){
		print UYT $_;
	};
	close(DFH);
	#system("cat make_fgts.sdf ring.sdf > a");
	#rename "a","make_fgts.sdf";
};
if(-e "fused_rings.sdf" && !-z "fused_rings.sdf"){
	open(DFH,"fused_rings.sdf");
	while(<DFH>){
		print UYT $_;
	};
	close(DFH);
	#system("cat make_fgts.sdf fused_rings.sdf > a");
	#rename "a","make_fgts.sdf";
};
if(-e "special.sdf" && !-z "special.sdf"){
	open(DFH,"special.sdf");
	while(<DFH>){
                print UYT $_;
        };
        close(DFH);
	#system("cat make_fgts.sdf special.sdf > a");
	#rename "a","make_fgts.sdf";
};
if(-e "linker.sdf" && !-z "linker.sdf"){
	open(DFH,"linker.sdf");
	while(<DFH>){
                print UYT $_;
        };
        close(DFH);
	#system("cat make_fgts.sdf linker.sdf > a");
	#rename "a","make_fgts.sdf";
};
if(-e "substituent.sdf" && !-z "substituent.sdf"){
	open(DFH,"substituent.sdf");
	while(<DFH>){
                print UYT $_;
        };
        close(DFH);
	#system("cat make_fgts.sdf substituent.sdf > a");
	#rename "a","make_fgts.sdf";
};
close(UYT);

if ($flagmend){
	print "unusual format of the sdf file : no \'\$\$\$\$\' \n"; 
	&sdf;
};


##############################################################################################


sub sdf{

	@copylign='';
	@copylign=@lignecor;

########## FRAGMENTS une Liaison

#	open(OUTC,">>fgt.sdf");
	foreach $k (1..$nbatom){
		if ($type[$k] ne 'H'){
			@get = split(' ',$bond[$k]);
			$long=@get/2;
#			printf OUTC "\n";
#			printf OUTC "  -ISIS-  08029913582D\n";
#			printf OUTC "\n";
			$longa=$long+1;
#			printf OUTC "%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$longa,$long;
#			printf OUTC "$lignecor[$k]";
			$m=0;
			foreach $l (1..$long){
#				printf OUTC "$lignecor[$get[$m]]";
				$m=$m+2;
			};
			$m=0;
			$y=1;
			foreach $l (1..$long){
				$x=1;
				$y=$y+1;
#				printf OUTC "%3s%3s%3s  0  0  0  0\n",$x,$y,$get[$m+1];
				$m=$m+2;
			};
#			printf OUTC "M  END\n";
#			printf OUTC "\n";
#			printf OUTC "\$\$\$\$\n";
		};
	};
#	close(OUTC);

	
########## FRAGMENTS (reconnaissance de cycles etc ....)

			print "Warning \'X\' dummy atom detected !\n" if($atom[$ig] eq 'X');
	$istratom=$nbatom;
	&cyclesdf;

#print "$nbc cycles\n";

#print "@copylign\n";
	foreach $lc (0..@cyclef-1){
		@copylig=split(' ',$cyclef[$lc]);
		foreach $copyl (0..@copylig-1){
			$copylign[$copylig[$copyl]]='';
		};
	};
#print "@copylign\n";

		
########## CYCLES JOINTS

		$nbcf=0;
		@cycleff='';
		@cyclefb='';

		foreach $lc (0..$nbc-2){
			foreach $lc2 (($lc+1)..$nbc-1){
				@get1= split(' ',$cyclef[$lc]);
				@get2= split(' ',$cyclef[$lc2]);
					foreach $lc3 (0..@get1-2){
						$couple1=$get1[$lc3].$blanc.$get1[$lc3+1];
						$couple1b=$get1[$lc3+1].$blanc.$get1[$lc3];
						foreach $lc4 (0..@get2-2){
							$couple2=$get2[$lc4].$blanc.$get2[$lc4+1];
							if (($couple1 eq $couple2)||($couple1b eq $couple2)){
#print "bond commun $couple1 \n";
								$cycleff[$nbcf]=$lc.$blanc.$lc2 if ($cycleff[$nbcf] eq '');
								$cyclefb[$nbcf]=$cyclefb[$nbcf].$blanc.$couple1;
							}
							elsif (($get1[$lc3] eq $get2[$lc4])||($get1[$lc3] eq $get2[$lc4+1])){
								$cycleff[$nbcf]=$lc.$blanc.$lc2 if ($cycleff[$nbcf] eq '');
							}
							elsif (($get1[$lc3+1] eq $get2[$lc4])||($get1[$lc3+1] eq $get2[$lc4+1])){
								$cycleff[$nbcf]=$lc.$blanc.$lc2 if ($cycleff[$nbcf] eq '');
							};
						};
						$couple2=$get2[0].$blanc.$get2[@get2-1];
						if (($couple1 eq $couple2)||($couple1b eq $couple2)){
#print "bond commun $couple1 \n";
							$cycleff[$nbcf]=$lc.$blanc.$lc2 if ($cycleff[$nbcf] eq '');
							$cyclefb[$nbcf]=$cyclefb[$nbcf].$blanc.$couple1;
						};
					};
					$couple1=$get1[@get1-1].$blanc.$get1[0];
					$couple1b=$get1[0].$blanc.$get1[@get1-1];
					foreach $lc4 (0..@get2-2){
						$couple2=$get2[$lc4].$blanc.$get2[$lc4+1];
						if (($couple1 eq $couple2)||($couple1b eq $couple2)){
#print "bond commun $couple1 \n";
							$cycleff[$nbcf]=$lc.$blanc.$lc2 if ($cycleff[$nbcf] eq '');
							$cyclefb[$nbcf]=$cyclefb[$nbcf].$blanc.$couple1;
						};
					};
					$couple2=$get2[0].$blanc.$get2[@get2-1];
					if (($couple1 eq $couple2)||($couple1b eq $couple2)){
#print "bond commun $couple1 \n";
						$cycleff[$nbcf]=$lc.$blanc.$lc2 if ($cycleff[$nbcf] eq '');
						$cyclefb[$nbcf]=$cyclefb[$nbcf].$blanc.$couple1;
					};
				 $nbcf++ if ($cycleff[$nbcf] ne '');
			};
		};

#########

@cyclefp=@cyclef;
if ($nbcf > 0){
		foreach $lc2 (0..$nbcf-1){
#			print "cycleff $cycleff[$lc2]\n";
			@get = split(' ',$cycleff[$lc2]);
			$cyclefp[$get[0]]='';
			$cyclefp[$get[1]]='';
		};

		foreach $lc (0..$nbcf-2){
			if ($cycleff[$lc] ne ''){

			@get = split(' ',$cycleff[$lc]);
			$lc4=0;
			while ($lc4 <= @get-1){
#print "search $lc4\n";
			foreach $lc2 (($lc+1)..$nbcf-1){
				if ($cycleff[$lc2] ne ''){
				@get2 = split(' ',$cycleff[$lc2]);
				if ($get[$lc4] eq $get2[0]){
					$test=$blanc.$get2[1].$blanc;
					$testb=$blanc.$cycleff[$lc].$blanc;
					if ($testb=~/($test)/){
						$rien = 0;
						$cycleff[$lc2]='';
					}
					else{
						$cycleff[$lc]=$cycleff[$lc].$blanc.$get2[1];
						$cycleff[$lc2]='';
					};
				}
				elsif($get[$lc4] eq $get2[1]){
					$test=$blanc.$get2[0].$blanc;
					$testb=$blanc.$cycleff[$lc].$blanc;
					if ($testb=~/($test)/){
						$rien = 0;
						$cycleff[$lc2]='';
					}
					else{
						$cycleff[$lc]=$cycleff[$lc].$blanc.$get2[0];
						$cycleff[$lc2]='';
					};
				};
				};
			};
#print "$cycleff[$lc]\n";
			$lc4++;
			@get = split(' ',$cycleff[$lc]);
			};
		};
		};

#		foreach $lc2 (0..$nbcf-1){
#			print "apres fusion $cycleff[$lc2]\n";
#		};

		foreach $lc (0..$nbcf-2){
			foreach $lc2 (($lc+1)..$nbcf-1){
				@get = split(' ',$cycleff[$lc]);
				@get2 = split(' ',$cycleff[$lc2]);
				$nbzero=0;
				foreach $lc3 (0..@get){
					foreach $lc4 (0..@get2){
						if (@get == @get2){
							if ($get2[$lc4] ne ''){
								if ($get[$lc3] eq $get2[$lc4]){
									$get2[$lc4]='';
									$nbzero ++;
								};
							};
						};
					};
				};
				$cycleff[$lc2]= '' if ($nbzero == @get2);
			};
		};


		$nbv=0;
		@reduit='';
		foreach $lc2 (0..$nbcf-1){
			$reduit[$nbv]=$cycleff[$lc2] if ($cycleff[$lc2] ne '');
			$nbv++ if ($cycleff[$lc2] ne '')
		};
		@cycleff='';
		@cycleff=@reduit;
		$nbcf=@cycleff;
#		foreach $lc2 (0..$nbcf-1){
#			print "$cycleff[$lc2]\n";
#		};


############# ELIMINATION DOUBLES

		$fg=0;
		@cyclefg='';
		@cyclefg2='';
		foreach $lc (0..$nbcf-1){
			@get = split(' ',$cycleff[$lc]);
			$longget=@get;

			@get1= split(' ',$cyclef[$get[0]]);
			$long1=@get1;

			@put='';
			$cyclefg2[$fg]=$cycleff[$lc];
			foreach $rg (0..@get1-1){
				$cyclefg[$fg]=$cyclefg[$fg].$blanc.$get1[$rg];
				$put[$rg]=$get1[$rg];
			};
			foreach $il (1..$longget){
				@get2= split(' ',$cyclef[$get[$il]]);
				$long2=@get2;

				foreach $rg (0..@put-1){
				foreach $rg2 (0..@get2-1){
					$get2[$rg2]='' if ($put[$rg] == $get2[$rg2]);
				};
				};

				$rh=@put;
				foreach $rg (0..@get2-1){
					if ($get2[$rg] ne ''){
						$put[$rh]=$get2[$rg];
						$cyclefg[$fg]=$cyclefg[$fg].$blanc.$get2[$rg];
						$rh++;
					};
				};
			};
			$fg ++
		};

		foreach $lc (0..$fg-2){
			foreach $lc2 (($lc+1)..$fg-1){
				@get = split(' ',$cyclefg[$lc]);
				@get2 = split(' ',$cyclefg[$lc2]);
				$nbzero=0;
				foreach $lc3 (0..@get){
					foreach $lc4 (0..@get2){
						if (@get == @get2){
							if ($get2[$lc4] ne ''){
								if ($get[$lc3] eq $get2[$lc4]){
									$get2[$lc4]='';
									$nbzero ++;
								};
							};
						};
					};
				};
				$cyclefg[$lc2]= '' if ($nbzero == @get2);
			};
		};

		$nbv=0;
		@reduit='';

		foreach $lc2 (0..$fg-1){
			$reduit[$nbv]=$cyclefg2[$lc2] if ($cyclefg[$lc2] ne '');
			$nbv++ if ($cyclefg[$lc2] ne '');
		};

		@cycleff='';
		@cycleff=@reduit;
		$nbcf=@cycleff;

#		foreach $lc2 (0..$nbcf-1){
#			print"$cycleff[$lc2]\n";
#		};
#	print "$nbcf systemes fusionnes\n";
#	print "\n";

		
################################################################		
# SPECIAL CASES: ring systems linked by double bonds
# eg: drug CYPROHEPTADINE	
		
$ligneb_special="";
foreach $lc (0..$nbcf-1){

	
	$pointcycleother='';
	foreach $lc2 (0..$nbc-1){
		if ($cyclefp[$lc2] ne ''){
			$pointcycleother=$pointcycleother." $cyclefp[$lc2] ";
		};
	};
	#print "other atoms cycles simple $pointcycleother\n";
			
	$pointcycle="";
	@get=split(' ',$cycleff[$lc]);
	foreach $lcp (0..@get-1){
		$pointcycle=$pointcycle." ".$cyclef[$get[$lcp]]." ";
	};	
	#print "Actuel cycles fusionnes $cycleff[$lc] => -$pointcycle-\n";
	
	$pointcycleotherf='';
	foreach $lc2 (0..$nbcf-1){
		if($lc != $lc2){
			@getp=split(' ',$cycleff[$lc2]);
			foreach $lcp (0..@getp-1){
				$pointcycleotherf=$pointcycleotherf." ".$cyclef[$getp[$lcp]]." ";
			};	
		};
	};
	#print "other atoms cycles fusionnes $pointcycleotherf\n";

	@get=split(' ',$pointcycle);
	$lcp=0;
	while($lcp < @get){
		$doublerings=0;
		@get2b = split(' ',$listb[$get[$lcp]]);
		foreach $lcp2 (0..@get2b-1){
			
			if($pointcycleotherf=~/ $get2b[$lcp2] /){
				@get3b = split(' ',$typeb[$get[$lcp]]);
				if($get3b[$lcp2]==2){
					#print "atom $get[$lcp] double bond with another fused ring system (atom $get2b[$lcp2])\n";
				foreach $lc2 (0..$nbcf-1){
					if($lc2 != $lc){
						@get4b = split(' ',$cycleff[$lc2]);
						foreach $lc3 (0..@get4b-1){
							if($cyclef[$get4b[$lc3]]=~/ $get2b[$lcp2] / && $cycleff[$lc]!~/ $get4b[$lc3] /){
								$cycleff[$lc]=$cycleff[$lc]." ".$cycleff[$lc2];
								$cycleff[$lc2]="";
								#print "add contents $cycleff[$lc2] => $cycleff[$lc]\n";
								$doublerings=1;
								$ligneb_special=$ligneb_special." $get[$lcp]-$get2b[$lcp2] ";
								#print "$ligneb_special\n";
							};
						};
					};
				};
				};
			};
			
			if($pointcycleother=~/ $get2b[$lcp2] /){
				@get3b = split(' ',$typeb[$get[$lcp]]);
				if($get3b[$lcp2]==2){
					#print "atom $get[$lcp] double bond with another single ring system (atom $get2b[$lcp2])\n";
				foreach $lc2 (0..$nbc-1){
					if($cyclef[$lc2]=~/ $get2b[$lcp2] /){
						$cycleff[$lc]=$cycleff[$lc]." ".$lc2;
						$cyclefp[$lc2]="";
						#print "add $lc2 => $cycleff[$lc]\n";
						$doublerings=1;
						$ligneb_special=$ligneb_special." $get[$lcp]-$get2b[$lcp2] ";
						#print "$ligneb_special\n";
					};	
				};
				};
			};
			
		};
		if($doublerings){
			$pointcycleother='';
			foreach $lc2 (0..$nbc-1){
				if ($cyclefp[$lc2] ne ''){
					$pointcycleother=$pointcycleother." $cyclefp[$lc2] ";
				};
			};
			#print "other atoms cycles simple $pointcycleother\n";
			
			$pointcycle="";
			@get=split(' ',$cycleff[$lc]);
			foreach $lcp (0..@get-1){
				$pointcycle=$pointcycle." ".$cyclef[$get[$lcp]]." ";
			};
			#print "Actuel cycles fusionnes $cycleff[$lc] => -$pointcycle-\n";
			
			$pointcycleotherf='';
			foreach $lc2 (0..$nbcf-1){
				if($lc != $lc2){
					@getp=split(' ',$cycleff[$lc2]);
					foreach $lcp (0..@getp-1){
						$pointcycleotherf=$pointcycleotherf." ".$cyclef[$getp[$lcp]]." ";
					};
				};
			};
			#print "other atoms cycles fusionnes $pointcycleotherf\n";
			
			@get=split(' ',$pointcycle);	
		};
		$lcp++;	
	};	
	
};	

		
############# SDF DES CYCLES FUSIONNE
#
	open(OUTD,">>special_tmp.sdf");
	open(OUTC,">>cyclefus_tmp.sdf");

	$nbatomar=0;
	$nbringar=0;

	foreach $lc (0..$nbcf-1){
	if($cycleff[$lc] ne ""){	
	#print "systeme no $lc\n";
	                $pointcycle='';
	                @new_coord='';
			$point='';

			@get = split(' ',$cycleff[$lc]);
			@getp = split(' ',$cyclefb[$lc]);
			

#######################################################
# reordonne les cycles du plus petit au plus grand sinon peut generer un manque de liaison
			@getcopy=@get;
			@getsort='';
			@lgetcyclef='';
			foreach $lcp (0..@get-1){
				@cycleftemp=split(' ',$cyclef[$get[$lcp]]);
				$longcycleftemp=@cycleftemp;

				if($longcycleftemp==5 || $longcycleftemp==6){	
				#print "ring @cycleftemp\n";
				$nbatomarpartial=0;	
				foreach $nbringi (0..@cycleftemp-1){
					$narp=" ".$typeb[$cycleftemp[$nbringi]]." ";
					$nocc1=$narp=~s/1/1/g;
					$nocc2=$narp=~s/2/2/g;
					if($type[$cycleftemp[$nbringi]] eq 'C'){
						$nbatomarpartial++ if($nocc1==2 && $nocc2==1);
					}
					elsif($type[$cycleftemp[$nbringi]] eq 'O'){
						$nbatomarpartial++ if($nocc1==2 && $longcycleftemp==5);
					}
					elsif($type[$cycleftemp[$nbringi]] eq 'N'){
						if($nocc1==1 && $nocc2==1){
							$nbatomarpartial++;
						}
						elsif($nocc1==3 && $longcycleftemp==5){
							$nbatomarpartial++;
						};
					}	
					elsif($type[$cycleftemp[$nbringi]] eq 'S'){
						$nbatomarpartial++ if($nocc1==2 && $longcycleftemp==5);
					};	
					if($nbringi == @cycleftemp-1){
						$nbringar++ if($nbatomarpartial == $longcycleftemp);
					};		
				};
				#print "\tring$nbringar\n";
				};
				
				$lgetcyclef[$lcp]=$longcycleftemp;
			};
#print "longue @lgetcyclef\n";


			$leplusptt=$lgetcyclef[0];
			$leplusptti=0;
			$getsorti=0;
			foreach $lcp2 (0..@get-1){
				foreach $lcp (0..@get-1){
					if($leplusptt > $lgetcyclef[$lcp] && $getcopy[$lcp] ne ''){
						$leplusptt=$lgetcyclef[$lcp];
						$leplusptti=$lcp;
					};
				};
				$getsort[$getsorti]=$getcopy[$leplusptti];
				$getcopy[$leplusptti]='';
				$getsorti++;

				foreach $lcp (0..@get-1){
					if($getcopy[$lcp] ne ''){
						$leplusptt=$lgetcyclef[$lcp];
						$leplusptti=$lcp;
					};
				};
			};

			@get=@getsort;


#print "@getsort\n";

#######################################################

#print "couples/liaisons @getp\n";
			$longget=@get;
			$long=@getp/2; #@getp n'est plus utilise par la suite
			
			foreach $lcp (0..@get-1){
				$pointcycle=$pointcycle." $cyclef[$get[$lcp]] ";
			};
			@addh='';
			@addh_bond='';
			$addhi=0;
			$listatom=$pointcycle;
			#print "atomes $pointcycle\n";
			#print "Fusion de $longget cycles $cycleff[$lc] ($long liaisons communes) \n";
## cycle 1
			@get1= split(' ',$cyclef[$get[0]]);
			#print "cycle1 $cyclef[$get[0]]\n";
			$long1=@get1;

			@put='';
			@ligneb='';
			$lb=0;

			@put2='';
			$p2=0;
			@ligneb2='';
			$lb2=0;

			@put3='';
			$p3=0;
			@ligneb3='';
			$lb3=0;

			foreach $rg (0..@get1-1){

				$put[$rg]=$get1[$rg];
				#print "\natom cycle $get1[$rg] [put = @put] [put2 = @put2] [put3 = @put3]\n";

				$situ=10000; #par defaut on n'ajoute pas l'atome

				@getb = split(' ',$listb[$put[$rg]]);
				#print "put $listb[$put[$rg]]\n";
				
				
				$narp=" ".$typeb[$put[$rg]]." ";
				$nocc1=$narp=~s/1/1/g;
				$nocc2=$narp=~s/2/2/g;
				if($type[$put[$rg]] eq 'C'){
					$nbatomar++ if($nocc1==2 && $nocc2==1);
                                }
                                elsif($type[$put[$rg]] eq 'O'){
					$nbatomar++ if($nocc1==2 && $long1==5);
		                }
				elsif($type[$put[$rg]] eq 'N'){
					if($nocc1==1 && $nocc2==1){
						$nbatomar++;
					}
					elsif($nocc1==3 && $long1==5){
						$nbatomar++;
					};
				}
				elsif($type[$put[$rg]] eq 'S'){
					$nbatomar++ if($nocc1==2 && $long1==5);
				};		
				
				foreach $rfb (0..@getb-1){
			
					$situ=$rfb if ($get1[$rg+1] == $getb[$rfb]); # si c'est un voisin accroche
					@getb2 = split(' ',$listb[$getb[$rfb]]);
					#print "getb atom $getb[$rfb] = $listb[$getb[$rfb]]\n";

					if (@getb2 == 1){ # si atome voisin a une seule liaison ex : hydrogene on le met dans put2
						$put2[$p2]=$getb[$rfb];
						@get3b = split(' ',$typeb[$put[$rg]]);
						$x=$rg+1;
						$y=$p2;
						$z="  0  0  0  0";
						$ligneb2[$lb2]=sprintf"%3s%3s%3s  0  0  0  0",$x,$y,$get3b[$rfb];
#print "1 b2 $ligneb2[$lb2]\n";
						$lb2++;
						$p2++;
					}
					else{
						$flagb=1;
						foreach $b2 (0..@getb2-1){
							@getb3 = split(' ',$listb[$getb2[$b2]]);
							if ($getb2[$b2] != $put[$rg]){
								$flagb=0 if(@getb3 > 1);
							};
							#print "n'ajoute pas $getb2[$b2] car c'est un fragment (a plusieurs substituants) $put[$rg] \n" if($flagb == 0);
						};
						if ($flagb){ # c'est pas exemple un C=N-H avec le C dans le cycle fusionne
							$put2[$p2]=$getb[$rfb];
							@get3b = split(' ',$typeb[$put[$rg]]);
							$x=$rg+1;
							$y=$p2;
							$z="  0  0  0  0";
							$ligneb2[$lb2]=sprintf"%3s%3s%3s  0  0  0  0",$x,$y,$get3b[$rfb];
#print "b2 $ligneb2[$lb2]\n";
							$lb2++;
							$p2++;

							foreach $b2 (0..@getb2-1){
								if ($getb2[$b2] != $put[$rg]){
									$put3[$p3]=$getb2[$b2];
#print "$put3[$p3]\n";
									@get4b = split(' ',$typeb[$getb[$rfb]]);
									$x=$p2-1;
									$y=$p3;
									$z="  0  0  0  0";
									$ligneb3[$lb3]=sprintf"%3s%3s%3s  0  0  0  0",$x,$y,$get4b[$b2];
#print "b3 $ligneb3[$lb3]\n";
									$lb3++;
									$p3++;
								};
							};
						};
					};
				};
				if ($situ != 10000){
					@get3b = split(' ',$typeb[$put[$rg]]);
					$x=$rg+1;
					$y=$x+1;
					$z="  0  0  0  0";
					$ligneb[$lb]=sprintf"%3s%3s%3s  0  0  0  0",$x,$y,$get3b[$situ];
					#print "ici $ligneb[$lb]\n";
#print "b4 $ligneb[$lb]\n";
					$lb++;
				};

			};
			$situ=10000;
			@getb = split(' ',$listb[$put[0]]);
			foreach $rfb (0..@getb-1){
				$situ=$rfb if ($get1[@get1-1] == $getb[$rfb]);
			};
			if ($situ != 10000){
				@get3b = split(' ',$typeb[$put[0]]);
				$x=1;
				$y=@get1;
				$z="  0  0  0  0";
				$ligneb[$lb]=sprintf"%3s%3s%3s  0  0  0  0",$x,$y,$get3b[$situ];
				#print "ici2 $ligneb[$lb]\n";
				#print "$ligneb[$lb]\n";
				$lb++;
			};
#print "@put\n";
			$nbfusion=1;
			foreach $il (1..($longget-1)){
				$putanc=@put;
				&fusion($get[$il]);
				$nbfusion++ if ($putanc < @put);
			};

##################################################### PRINT cyclefus.sdf


#			print "put @put\n";
			$longa=$rh;
			$longb=$lb;
			#print "$longa atomes et $longb bond\n";
			$nhetat=0;
			$nar=0;

			#printf OUTC "\n";
			#printf OUTC "  -ISIS-  08029913582D\n";
			#printf OUTC "\n";
			#printf OUTC "%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$longa,$longb;

			$new_coord[0]="\n";
			$new_coord[1]="  -ISIS-  08029913582D\n";
			$new_coord[2]="\n";
			$new_coord[3]=sprintf"%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$longa,$longb;
			$longa2=$longa;
			$longb2=$longb;

			@correspondance="";
			
			$rgx=0;
			$rgy=0;
			$rgz=0;
			$rgnb=0;
			foreach $l (0..$rh-1){
				#printf OUTC "$lignecor[$put[$l]]";
				#print "$lignecor[$put[$l]]";
				$new_coord[4]=$new_coord[4]."$lignecor[$put[$l]]";
				#print "coord $l -> old $put[$l]\n";
				$correspondance[$put[$l]]=$l+1;

				$nhetat++ if (($lignecor[$put[$l]]=~/O/)||($lignecor[$put[$l]]=~/N/)||($lignecor[$put[$l]]=~/S/));
				$differentielx++ if($lignecor[$put[$l]]=~/ X /);

				@get7=split(' ',$lignecor[$put[$l]]);
				# Geometric center
				if($lignecor[$put[$l]]!~/ H / && $lignecor[$put[$l]]!~/ X /){
					$rgx=$rgx+$get7[0];
					$rgy=$rgy+$get7[1];
					$rgz=$rgz+$get7[2];
					$rgnb++;
				};

				$flag_point=0;
				@get2pt = split(' ',$listb[$put[$l]]);
				@get2pttype = split(' ',$typeb[$put[$l]]);

				foreach $lp (0..@get2pt-1){
					$point_atome=$get2pt[$lp];
					$get2pt[$lp]=" ".$get2pt[$lp]." ";

				 	if($pointcycle!~/$get2pt[$lp]/){

				 			$lpoint=$l+1;
				 			@get3 = split(' ',$lignecor[$point_atome]);
				 		
				 			if($type[$put[$l]] eq "P" && $type[$point_atome] eq "O"){
				 				
				 				@get4=split(' ',$listb[$put[$l]]);
				 				@get5=split(' ',$typeb[$put[$l]]);
				 				foreach $ddp (0..@get4-1){
				 			        	$ddpref=$ddp if($get4[$ddp]==$point_atome);
				 				};
				 			
				 				$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get3[0],$get3[1],$get3[2],$blanc,$get3[3],$blanc,$get3[4],$get3[5],$get3[6];
				 			
								# Geometric center
								if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
			                                        	$rgx=$rgx+$get3[0];
                        			                	$rgy=$rgy+$get3[1];
 	                                  				$rgz=$rgz+$get3[2];
								$rgnb++;
								};

				 				$y=$longa+$addhi+1;
				 				$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$lpoint,$y,$get5[$ddpref];
				 			
				 				$addhi++;
				 				$longa2++;
				 				$longb2++;
				 				
				 				$copylign[$point_atome]='';
				 				$listatom=$listatom." ".$point_atome." ";
				 				
				 				if($get5[$ddpref]==1){
				 					@get4=split(' ',$listb[$point_atome]);
				 					$ddpref='';
				 					foreach $ddp (0..@get4-1){
				 			        		$ddpref=$ddp if($type[$get4[$ddp]] ne 'P');
				 					};
				 					if($ddpref ne ''){
				 						@get6 = split(' ',$lignecor[$get4[$ddpref]]);
				 						
				 						if($get6[3] ne "H"){
				 							if($point eq ''){
												$point=$y;
												$pointo=$point_atome."_".$get4[$ddpref];
											}
											else{
												$point=$point."-".$y;
												$pointo=$pointo."-".$point_atome."_".$get4[$ddpref];
											};
										
				 							$get6[3]="X"; #si on veux un dummy atom
				 						}
				 						else{
				 							$copylign[$get4[$ddpref]]='';		
				 						};
				 						
				 						$y2=$longa+$addhi+1;
				 						$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get6[0],$get6[1],$get6[2],$blanc,$get6[3],$blanc,$get6[4],$get6[5],$get6[6];
				 			                       
										# Geometric center
										if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
											$rgx=$rgx+$get6[0];
											$rgy=$rgy+$get6[1];
											$rgz=$rgz+$get6[2];
											$rgnb++;
										};

									 	$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$y,$y2,1;
				 			                        $addhi++;
				 						$longa2++;
				 						$longb2++;
				 						
				 						$copylign[$get4[$ddpref]]='';
				 					};
				 				};
				 			}
				 			elsif($get2pttype[$lp]==2){ #recherche si double liaison
				 				
				 				@get4=split(' ',$listb[$put[$l]]);
				 				@get5=split(' ',$typeb[$put[$l]]);
				 				foreach $ddp (0..@get4-1){
				 			        	$ddpref=$ddp if($get4[$ddp]==$point_atome);
				 				};
				 			
				 				$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get3[0],$get3[1],$get3[2],$blanc,$get3[3],$blanc,$get3[4],$get3[5],$get3[6];
				 			
								# Geometric center
								if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
									$rgx=$rgx+$get3[0];
									$rgy=$rgy+$get3[1];
									$rgz=$rgz+$get3[2];
									$rgnb++;
								};

				 				$y=$longa+$addhi+1;
				 				$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$lpoint,$y,$get5[$ddpref];
				 			        $voisin2bis=$y;
				 				$addhi++;
				 				$longa2++;
				 				$longb2++;
				 				
				 				$copylign[$point_atome]='';
				 			        $listatom=$listatom." ".$point_atome." ";
				 			
								$addh_voisin2=" $point_atome ";
								$addh_voisin2i=" $y ";
								$voisin2=$point_atome;
								$addhi_old=$addhi;
								$continuedouble=1;
								while($continuedouble){
		
									@get4=split(' ',$listb[$voisin2]);
									@get5=split(' ',$typeb[$voisin2]);
									foreach $ddp (0..@get5-1){
										if($get5[$ddp]>1 && $listatom!~/ $get4[$ddp] /){
			
											@get6 = split(' ',$lignecor[$get4[$ddp]]);	 			 		
											$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get6[0],$get6[1],$get6[2],$blanc,$get6[3],$blanc,$get6[4],$get6[5],$get6[6];
				 			
											# Geometric center
											if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
												$rgx=$rgx+$get6[0];
												$rgy=$rgy+$get6[1];
												$rgz=$rgz+$get6[2];
												$rgnb++;
											};

											$y=$longa+$addhi+1;
											$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$voisin2bis,$y,$get5[$ddp];
				 			
											$addhi++;
											$longa2++;
											$longb2++;
				
											$addh_voisin2=$addh_voisin2.$get4[$ddp]." ";
											$addh_voisin2i=$addh_voisin2i.$y." ";
											$copylign[$get4[$ddp]]='';
											$listatom=$listatom." ".$get4[$ddp]." ";			 		
										};
									};
									if($addhi>$addhi_old){     # suppose (et c'est toujours vrai) que l'atome n'a pas plus de 2 doubles liaisons
										$addhi_old=$addhi;
										$voisin2=$get4[$ddp];
										$voisin2bis=$y;
									}
									else{
										$continuedouble=0;
									};
		
								};
	
								#Les H ou X a ajouter
								@get6=split(' ',$addh_voisin2);	
								@get8=split(' ',$addh_voisin2i);			 				
								foreach $hil (0..@get6-1){
									@get9=split(' ',$listb[$get6[$hil]]);
									foreach $ikhp (0..@get9-1){
	       									@get7 = split(' ',$lignecor[$get9[$ikhp]]);
	       								        if($listatom!~/ $get9[$ikhp] /){
	       										if($get7[3] ne "H"){
				 								if($point eq ''){
													$point=$get8[$hil];
													$pointo=$get6[$hil]."_".$get9[$ikhp];
												}
												else{
													$point=$point."-".$get8[$hil];
													$pointo=$pointo."-".$get6[$hil]."_".$get9[$ikhp];
												};
										
				 								$get7[3]="X"; #si on veux un dummy atom
				 							}
				 							else{
				 								$copylign[$get9[$ikhp]]='';	
				 							};
				 				
				 			
				 							$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get7[0],$get7[1],$get7[2],$blanc,$get7[3],$blanc,$get7[4],$get7[5],$get7[6];
				 							
											# Geometric center
											if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
                                                                                        	$rgx=$rgx+$get7[0];
                                                                                        	$rgy=$rgy+$get7[1];
												$rgz=$rgz+$get7[2];
												$rgnb++;
											};

				 							$y=$longa+$addhi+1;
				 							$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$get8[$hil],$y,1;
				 			
				 							$addhi++;
				 							$longa2++;
				 							$longb2++;
	                                                                    };
	           							};
								};			 				
				 			
				 			}
				 			else{	
				 				if($get3[3] ne "H"){
				 					if($point eq ''){
										$point=$lpoint;
										$pointo=$put[$l]."_".$point_atome;
									}
									else{
										$point=$point."-".$lpoint;
										$pointo=$pointo."-".$put[$l]."_".$point_atome;
									};
										
				 					$get3[3]="X"; #si on veux un dummy atom
				 				}
				 				else{
				 					$copylign[$point_atome]='';
				 				};
				 				
				 			
				 				$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get3[0],$get3[1],$get3[2],$blanc,$get3[3],$blanc,$get3[4],$get3[5],$get3[6];
				 			
								# Geometric center
								if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
									$rgx=$rgx+$get3[0];
									$rgy=$rgy+$get3[1];
									$rgz=$rgz+$get3[2];
									$rgnb++;
								};

				 				$y=$longa+$addhi+1;
				 				$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$lpoint,$y,1;
				 			
				 				$addhi++;
				 				$longa2++;
				 				$longb2++;
				 			
				 			};
				 		
				 			
				 	};#if not in pointcycle
				 	
				};
			};

			$new_coordspecial="";	
			if ($ligneb_special ne ""){
				@getspecial=split(' ',$ligneb_special);
				foreach $lcp (0..@getspecial-1){
					@getspecial2=split('-',$getspecial[$lcp]);
					$aspecial="";
					$bspecial="";
					$aspecial=$correspondance[$getspecial2[0]];
					$bspecial=$correspondance[$getspecial2[1]];
					if($aspecial ne "" && $bspecial ne ""){
						#print "new $aspecial $bspecial\n";
						$new_coordspecial=$new_coordspecial.sprintf"%3s%3s%3s  0  0  0  0\n",$aspecial,$bspecial,2;
						$longb2++;
						$getspecial[$lcp]="";
					};	
				};
			};
			#print "new:\n $new_coordspecial\n";
			
			foreach $l (0..$lb-1){
				#printf OUTC "$ligneb[$l]\n";
				#printf "ligne $ligneb[$l]\n";
				$new_coord[5]=$new_coord[5]."$ligneb[$l]\n";
				#print "$ligneb[$l]\n";

				#print "cc $new_coord[5]\n";
				@getnar=split(' ',$ligneb[$l]);
				$nar++ if ($getnar[2] == 2);
			};
			 
			#printf OUTC "M  END\n";
			#printf OUTC ">  <ncycles> (1)\n";
			#printf OUTC "$nbfusion\n";
			#printf OUTC "\n";
			#printf OUTC ">  <natomescycle> (1)\n";
			#printf OUTC "$longa\n";
			#printf OUTC "\n";
			#printf OUTC ">  <cyclearom> (1)\n";
			#printf OUTC "$nar\n";
			#printf OUTC "\n";
			#printf OUTC ">  <nhetat> (1)\n";
			#printf OUTC "$nhetat\n";
			#printf OUTC "\n";
			#printf OUTC "\$\$\$\$\n";
		
		if($new_coordspecial eq ''){
			print OUTC "$new_coord[0]";
			print OUTC "$new_coord[1]";
			print OUTC "$new_coord[2]";
			$new_coord[3]=sprintf"%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$longa2,$longb2;
			$differentiel=$differentiel+$longb2;
			print OUTC "$new_coord[3]";
			print OUTC "$new_coord[4]";
			foreach $rf (0..@addh-1){
				print OUTC $addh[$rf];
				$differentielx++ if($addh[$rf]=~/ X /);
			};
			print OUTC "$new_coord[5]";
			print OUTC "$new_coordspecial" if($new_coordspecial ne "");
			
			#print  "$new_coord[5]";
			foreach $rf (0..@addh-1){
				print OUTC $addh_bond[$rf];
				#print $addh_bond[$rf];
			};
				
			if($rgnb != 0){
				$rgx=sprintf"%5.3f",$rgx/$rgnb;
				$rgy=sprintf"%5.3f",$rgy/$rgnb;
				$rgz=sprintf"%5.3f",$rgz/$rgnb;
			};

			print OUTC "M  END\n";
			print OUTC "> <ncycles>\n";
			print OUTC "$nbfusion\n";
			print OUTC "\n";
			print OUTC "> <natomescycle>\n";
			print OUTC "$longa\n";
			print OUTC "\n";
			print OUTC "> <cyclearom>\n";
			print OUTC "$nar\n";
			print OUTC "\n";
			print OUTC "> <nhetat>\n";
			print OUTC "$nhetat\n";
			print OUTC "\n";
			if($nbringar > 0){
				print OUTC "> <nbringar>\n";
				print OUTC "$nbringar\n";
				print OUTC "\n";
			};
			#if(($nbfusion == 2 && $longa == 10 && $nar == 5) || ($nbfusion == 2 && $longa == 9 && $nar ==4)){
			if($longa == $nbatomar){
				print OUTC "> <ARcenter>\n";
				print OUTC "$rgx $rgy $rgz\n";
				print OUTC "\n";
			}
			else{
				print OUTC "> <LIPcenter>\n";
				print OUTC "$rgx $rgy $rgz\n";
				print OUTC "\n";
			};

			if($point ne ''){
				print OUTC "> <POINTS>\n";
				print OUTC "$point\n";
				print OUTC "\n";
				print OUTC "> <ORIGIN>\n";
				print OUTC "$pointo\n";
				print OUTC "\n";
			};
			print OUTC "$lignedata" if($lignedata ne "");
			print OUTC "\$\$\$\$\n";
			
		}
		else{
			print OUTD "$new_coord[0]";
			print OUTD "$new_coord[1]";
			print OUTD "$new_coord[2]";
			$new_coord[3]=sprintf"%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$longa2,$longb2;
			$differentiel=$differentiel+$longb2;
			print OUTD "$new_coord[3]";
			print OUTD "$new_coord[4]";
			foreach $rf (0..@addh-1){
				print OUTD $addh[$rf];
				$differentielx++ if($addh[$rf]=~/ X /);
			};
			print OUTD "$new_coord[5]";
			print OUTD "$new_coordspecial" if($new_coordspecial ne "");
			
			foreach $rf (0..@addh-1){
				print OUTD $addh_bond[$rf];
			};	
			
			if($rgnb != 0){
				$rgx=sprintf"%5.3f",$rgx/$rgnb;
				$rgy=sprintf"%5.3f",$rgy/$rgnb;
				$rgz=sprintf"%5.3f",$rgz/$rgnb;
			};	
		
			print OUTD "M  END\n";
			print OUTD "> <ncycles>\n";
			print OUTD "$nbfusion\n";
			print OUTD "\n";
			print OUTD "> <natomescycle>\n";
			print OUTD "$longa\n";
			print OUTD "\n";
			print OUTD "> <cyclearom>\n";
			print OUTD "$nar\n";
			print OUTD "\n";
			print OUTD "> <nhetat>\n";
			print OUTD "$nhetat\n";
			print OUTD "\n";
			if($nbringar > 0){
				print OUTD "> <nbringar>\n";
				print OUTD "$nbringar\n";
				print OUTD "\n";
			};	

			#if(($nbfusion == 2 && $longa == 10 && $nar == 5) || ($nbfusion == 2 && $longa == 9 && $nar ==4)){
			if($nbatomar == $longa){
				print OUTD "> <ARcenter>\n";
				print OUTD "$rgx $rgy $rgz\n";
				print OUTD "\n";
			}
			else{
				print OUTD "> <LIPcenter>\n";
				print OUTD "$rgx $rgy $rgz\n";
				print OUTD "\n";
			};

			if($point ne ''){
				print OUTD "> <POINTS>\n";
				print OUTD "$point\n";
				print OUTD "\n";
				print OUTD "> <ORIGIN>\n";
				print OUTD "$pointo\n";
				print OUTD "\n";
			};
			print OUTD "$lignedata" if($lignedata ne "");
			print OUTD "\$\$\$\$\n";
			
		};	
	};# if cycleff ne ""
	};
	close(OUTC);
	close(OUTD);
};


################################################################
## SPECIAL CASES: ring systems linked by double bonds
#

$ligneb_special="";
foreach $lc2 (0..$nbc-1){
        if ($cyclefp[$lc2] ne ''){
		
		$pointcycleother='';
		foreach $lc ($lc2+1..$nbc-1){
			if ($cyclefp[$lc] ne ''){
				$pointcycleother=$pointcycleother." $cyclefp[$lc] ";
			};
		};
		#print "other atoms cycles simple $pointcycleother\n";
		
                @get=split(' ',$cyclef[$lc2]);
		$lcp=0;
                while($lcp < @get){
			$doublerings=0;
                        @get2b = split(' ',$listb[$get[$lcp]]);
                        foreach $lcp2 (0..@get2b-1){
                            	if($pointcycleother=~/ $get2b[$lcp2] /){
                            		@get3b = split(' ',$typeb[$get[$lcp]]);
					if($get3b[$lcp2]==2){
						#print "atom $get[$lcp] double bond with another single ring system (atom $get2b[$lcp2])\n";
				     		foreach $lc3 (0..$nbc-1){
					     		if($cyclef[$lc3]=~/ $get2b[$lcp2] /){
								$cyclef[$lc2]=$cyclef[$lc2]." * ".$cyclef[$lc3];
								$cyclefp[$lc3]="";
								#print "add ring $cyclef[$lc3] => $cyclef[$lc2]\n";
						     		$doublerings=1;
						     		$ligneb_special=$ligneb_special." $get[$lcp]-$get2b[$lcp2] ";
								#print "$ligneb_special\n";
							};
						};	
			   		};			
				};
                 	};
			if($doublerings){
				@get=split(' ',$cyclef[$lc2]);

				$pointcycleother='';
				foreach $lc ($lc2+1..$nbc-1){
					if ($cyclefp[$lc] ne ''){
						$pointcycleother=$pointcycleother." $cyclefp[$lc] ";
					};
				};
				#print "other atoms cycles simple $pointcycleother\n";	
			};	
			$lcp++;
           	};
		
       };
};

############################################################ CYCLE EN SDF

open(OUTC,">>cycle_tmp.sdf");
open(OUTD,">>special_tmp.sdf");


		foreach $lc (0..$nbc-1){
		
			$longa=0;
		        @new_coord='';
		        $point='';
			@addh='';
			@addh_bond='';
			$addhi=0;
			
			@correspondance='';
			
			if ($cyclefp[$lc] ne ''){

				#printf OUTC "\n";
				#printf OUTC "  -ISIS-  08029913582D\n";
				#printf OUTC "\n";
				#printf OUTC "%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$long,$long;
				
				
				$new_coord[0]="\n";
				$new_coord[1]="  -ISIS-  08029913582D\n";
				$new_coord[2]="\n";
				$new_coord[3]=sprintf"%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$long,$long;

				$nhetat=0;
				$nar=0;
				$nbatomar=0;
				$nbringar=0;
				
				$rgx=0;
				$rgy=0;
				$rgz=0;
				$rgnb=0;
				
				$pointcycle=" $cyclef[$lc] ";
				#print "-$pointcycle-\n";
				$pointcycle=~s/\*//g;
				$listatom=$pointcycle;
				
				
			# special case separate by star
			$speciali=0; # previous nb atom in ring 
			@get = split(' ',$listatom);
			$specialib=@get; # total nb atom rings, print addh numbers after it
			@getstar= split('\*',$cyclef[$lc]);
			$longgetstar=@getstar;
			
			foreach $lc3 (0..@getstar-1){
				@get = split(' ',$getstar[$lc3]);
				$long=@get;
				$longa=$longa+$long;

				$nbatomarpartial=0;
				foreach $l (0..$long-1){
					
					$new_coord[4]=$new_coord[4]."$lignecor[$get[$l]]";
					#printf OUTC "$lignecor[$get[$l]]";
					$correspondance[$get[$l]]=$l+1+$speciali;
					
					$nhetat++ if (($lignecor[$get[$l]]=~/O/)||($lignecor[$get[$l]]=~/S/)||($lignecor[$get[$l]]=~/N/));
					$differentielx++ if($lignecor[$get[$l]]=~/ X /);

					@get7 = split(' ',$lignecor[$get[$l]]);
					
					# Geometric center
					if($lignecor[$get[$l]]!~/ H / && $lignecor[$get[$l]]!~/ X /){
						$rgx=$rgx+$get7[0];
						$rgy=$rgy+$get7[1];
						$rgz=$rgz+$get7[2];
						$rgnb++;
					};
						
					$flag_point=0;
					@get2 = split(' ',$listb[$get[$l]]);
					@get2pttype = split(' ',$typeb[$get[$l]]);
					
					
					$narp=" ".$typeb[$get[$l]]." ";
					$nocc1=$narp=~s/1/1/g;
					$nocc2=$narp=~s/2/2/g;
					if($type[$get[$l]] eq 'C'){
						$nbatomar++ if($nocc1==2 && $nocc2==1);
						$nbatomarpartial++ if($nocc1==2 && $nocc2==1);
					}
					elsif($type[$get[$l]] eq 'O'){
						$nbatomar++ if($nocc1==2 && $long==5);
						$nbatomarpartial++ if($nocc1==2 && $long==5);
					}
					elsif($type[$get[$l]] eq 'N'){
						if($nocc1==1 && $nocc2==1){
							$nbatomar++;
							$nbatomarpartial++;
						}	
						elsif($nocc1==3 && $long==5){
							$nbatomar++;
							$nbatomarpartial++;
						};	
					}
					elsif($type[$get[$l]] eq 'S'){
						$nbatomar++ if($nocc1==2 && $long==5);
						$nbatomarpartial++ if($nocc1==2 && $long==5);
					};
					if($l == $long-1){
						$nbringar++ if($nbatomarpartial == $long);
					};	
					
					#print "nocc1 $nocc1 nocc2 $nocc2\n";
					#print "$specialib atoms in ring: $nbatomar ar\n";
					
					foreach $lp (0..@get2-1){
				        	$point_atome=$get2[$lp];
						$get2[$lp]=" ".$get2[$lp]." ";
				 		
				 		if($pointcycle!~/$get2[$lp]/){
							#print "ici point $get2[$lp]\n";
							
				 			$lpoint=$l+1+$speciali;
				 			
				 			@get3 = split(' ',$lignecor[$point_atome]);
				 		
				 			if($type[$get[$l]] eq "P" && $type[$point_atome] eq "O"){
				 			
				 			        @get4=split(' ',$listb[$get[$l]]);
				 			        @get5=split(' ',$typeb[$get[$l]]);
				 			        foreach $ddp (0..@get4-1){
				 			        	$ddpref=$ddp if($get4[$ddp]==$point_atome);
				 			        };
				 			
				 				$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get3[0],$get3[1],$get3[2],$blanc,$get3[3],$blanc,$get3[4],$get3[5],$get3[6];
				 		
								# Geometric center
								if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
                                                                	$rgx=$rgx+$get3[0];
                                                                	$rgy=$rgy+$get3[1];
                                                                	$rgz=$rgz+$get3[2];
									$rgnb++;
								};
		
				 				$y=$addhi+1+$specialib;
								$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$lpoint,$y,$get5[$ddpref];
								#print "addbond $lpoint - $y type = $get5[$ddpref]\n";
				 				$addhi++;
				 				$longa++;
				 				$copylign[$point_atome]='';
				 				$listatom=$listatom." ".$point_atome." ";
				 				
				 				if($get5[$ddpref]==1){
				 					@get4=split(' ',$listb[$point_atome]);
				 					$ddpref='';
				 					foreach $ddp (0..@get4-1){
				 			        		$ddpref=$ddp if($type[$get4[$ddp]] ne 'P');
				 					};
				 					if($ddpref ne ''){
				 						@get6 = split(' ',$lignecor[$get4[$ddpref]]);
				 						
				 						if($get6[3] ne "H"){
				 							if($point eq ''){
												$point=$y;
												$pointo=$point_atome."_".$get4[$ddpref];
											}
											else{
												$point=$point."-".$y;
												$pointo=$pointo."-".$point_atome."_".$get4[$ddpref];
											};
										
				 							$get6[3]="X"; #si on veux un dummy atom
				 						}
				 						else{
				 							$copylign[$get4[$ddpref]]='';	
				 						};
				 						
				 						$y2=$addhi+1+$specialib;
				 						$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get6[0],$get6[1],$get6[2],$blanc,$get6[3],$blanc,$get6[4],$get6[5],$get6[6];
				 			                     
										# Geometric center
										if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
		         	                                                       	$rgx=$rgx+$get6[0];
                			                                                $rgy=$rgy+$get6[1];
                                			                                $rgz=$rgz+$get6[2];
											$rgnb++;
										};

										$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$y,$y2,1;
										#print "addbond2 $y - $y2 type = 1\n";
				 			                        $addhi++;
				 						$longa++;
				 						$copylign[$get4[$ddpref]]='';
				 					};
				 				};
				 			}
				 			elsif($get2pttype[$lp]==2){ #recherche si double liaison
				 				
				 				@get4=split(' ',$listb[$get[$l]]);
				 				@get5=split(' ',$typeb[$get[$l]]);
				 				foreach $ddp (0..@get4-1){
				 			        	$ddpref=$ddp if($get4[$ddp]==$point_atome);
				 				};
				 			
				 				$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get3[0],$get3[1],$get3[2],$blanc,$get3[3],$blanc,$get3[4],$get3[5],$get3[6];
				 			
								# Geometric center
								if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
                                                                	$rgx=$rgx+$get3[0];
                                                                	$rgy=$rgy+$get3[1];
                                                                	$rgz=$rgz+$get3[2];
									$rgnb++;
								};

				 				$y=$addhi+1+$specialib;
								$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$lpoint,$y,$get5[$ddpref];
								#print "addbond3 $lpoint - $y type = $get5[$ddpref]\n";
				 			        $voisin2bis=$y;
				 				$addhi++;
				 				$longa++;
				 				
				 				
				 				$copylign[$point_atome]='';
				 			        $listatom=$listatom." ".$point_atome." ";
				 			
								$addh_voisin2=" $point_atome ";
								$addh_voisin2i=" $y ";
								$voisin2=$point_atome;
								$addhi_old=$addhi;
								$continuedouble=1;
								while($continuedouble){
		
									@get4=split(' ',$listb[$voisin2]);
									@get5=split(' ',$typeb[$voisin2]);
									foreach $ddp (0..@get5-1){
										if($get5[$ddp]>1 && $listatom!~/ $get4[$ddp] /){
			
											@get6 = split(' ',$lignecor[$get4[$ddp]]);	 			 		
											$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get6[0],$get6[1],$get6[2],$blanc,$get6[3],$blanc,$get6[4],$get6[5],$get6[6];


											# Geometric center
											if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
			                                                        	        $rgx=$rgx+$get6[0];
                        			                        	                $rgy=$rgy+$get6[1];
                                                			                	$rgz=$rgz+$get6[2];
												$rgnb++;
				 							};

											$y=$addhi+1+$specialib;
											$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$voisin2bis,$y,$get5[$ddp];
											#print "addbond4 $voisin2bis - $y type = $get5[$ddp]\n";
											$addhi++;
											$longa++;
											;
				
											$addh_voisin2=$addh_voisin2.$get4[$ddp]." ";
											$addh_voisin2i=$addh_voisin2i.$y." ";
											$copylign[$get4[$ddp]]='';
											$listatom=$listatom." ".$get4[$ddp]." ";			 		
										};
									};
									if($addhi>$addhi_old){     # suppose (et c'est toujours vrai) que l'atome n'a pas plus de 2 doubles liaisons
										$addhi_old=$addhi;
										$voisin2=$get4[$ddp];
										$voisin2bis=$y;
									}
									else{
										$continuedouble=0;
									};
		
								};
	
								#Les H ou X a ajouter
								@get6=split(' ',$addh_voisin2);	
								@get8=split(' ',$addh_voisin2i);			 				
								foreach $hil (0..@get6-1){
									@get9=split(' ',$listb[$get6[$hil]]);
									foreach $ikhp (0..@get9-1){
	       									@get7 = split(' ',$lignecor[$get9[$ikhp]]);
	       								        if($listatom!~/ $get9[$ikhp] /){
	       										if($get7[3] ne "H"){
				 								if($point eq ''){
													$point=$get8[$hil];
													$pointo=$get6[$hil]."_".$get9[$ikhp];
												}
												else{
													$point=$point."-".$get8[$hil];
													$pointo=$pointo."-".$get6[$hil]."_".$get9[$ikhp];
												};
										
				 								$get7[3]="X"; #si on veux un dummy atom
				 				                        }
				 							else{
				 								$copylign[$get9[$ikhp]]='';
				 							};
				 			
				 							$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get7[0],$get7[1],$get7[2],$blanc,$get7[3],$blanc,$get7[4],$get7[5],$get7[6];
				 			
											# Geometric center
											if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
			                                                                	$rgx=$rgx+$get7[0];
                        			                                        	$rgy=$rgy+$get7[1];
                                                			                	$rgz=$rgz+$get7[2];
												$rgnb++;
											};

				 							$y=$addhi+1+$specialib;
											$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$get8[$hil],$y,1;
											#print "addbond5 $get8[$hil] - $y type = 1\n";
				 							$addhi++;
				 							$longa++;
				 							
	                                                                    };
	           							};
								};			 				
				 			
				 			}
				 			else{	
				 				if($get3[3] ne "H"){
				 					if($point eq ''){
										$point=$lpoint;
										$pointo=$get[$l]."_".$point_atome;
									}
									else{
										$point=$point."-".$lpoint;
										$pointo=$pointo."-".$get[$l]."_".$point_atome;
									};
										
				 					$get3[3]="X"; #si on veux un dummy atom
				 				}
				 				else{
				 					$copylign[$point_atome]='';
				 				};
				 				
				 				$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get3[0],$get3[1],$get3[2],$blanc,$get3[3],$blanc,$get3[4],$get3[5],$get3[6];
				 			
				 				# Geometric center
								if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
									$rgx=$rgx+$get3[0];
									$rgy=$rgy+$get3[1];
									$rgz=$rgz+$get3[2];
									$rgnb++;
								};

								$y=$addhi+1+$specialib;
								$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$lpoint,$y,1;
								#print "addbond6 $lpoint - $y type = 1\n";
								
				 				$addhi++;
				 				$longa++;
				 			};	
				 		
				 	
				 		};
				 	
					};
					
				};
#print "neww_coord\n$new_coord[4]\n";
#print "adh\n@addh\n";

				foreach $l (0..$long-2){
					#print "correspond $get[$l] => $correspondance[$get[$l]]\n";
					
					$situ=10000;
					$x=$l+1+$speciali;
					$y=$x+1;
					@get2 = split(' ',$listb[$get[$l]]);
					foreach $rf (0..@get2-1){
						$situ=$rf if ($get[$l+1] == $get2[$rf]);
					};
					if ($situ != 10000){
						@get3 = split(' ',$typeb[$get[$l]]);
						$new_coord[5]=sprintf"$new_coord[5]%3s%3s%3s  0  0  0  0\n",$x,$y,$get3[$situ];

						#print "bond $x - $y type = $get3[$situ]\n";
						
						#printf OUTC "%3s%3s%3s  0  0  0  0\n",$x,$y,$get3[$situ];
						$nar++ if ($get3[$situ] == 2);
					};
				};
				#print "correspond $get[$long-1] => $correspondance[$get[$long-1]]\n";

				$x=$long+$speciali;
				$y=1+$speciali;
				$situ=10000;
				@get2 = split(' ',$listb[$get[0]]);
				foreach $rf (0..@get2-1){
					$situ=$rf if ($get[$long-1] == $get2[$rf]);
				};
				if ($situ != 10000){
					@get3 = split(' ',$typeb[$get[0]]);
					$new_coord[5]=sprintf"$new_coord[5]%3s%3s%3s  0  0  0  0\n",$x,$y,$get3[$situ];

					#print "bond $x - $y type = $get3[$situ]\n";

					#printf OUTC "%3s%3s%3s  0  0  0  0\n",$x,$y,$get3[$situ];
					$nar++ if ($get3[$situ] == 2);
				};

				$speciali=$speciali+$long;
				#print "speciali= $speciali in neww_coord\n";
				#print "specialib= $specialib in addh\n";
				#print "total = longa = $longa atoms\n";

			};#end foreach separate by star


				$longb2=$longa;
				$new_coordspecial="";
				if ($ligneb_special ne ""){
					@getspecial=split(' ',$ligneb_special);
					foreach $lcp (0..@getspecial-1){
						@getspecial2=split('-',$getspecial[$lcp]);
						$aspecial="";
						$bspecial="";
						$aspecial=$correspondance[$getspecial2[0]];
						$bspecial=$correspondance[$getspecial2[1]];
						if($aspecial ne "" && $bspecial ne ""){
							#print "new $aspecial $bspecial\n";
							$new_coordspecial=$new_coordspecial.sprintf"%3s%3s%3s  0  0  0  0\n",$aspecial,$bspecial,2;
							$longb2++;
							$getspecial[$lcp]="";
						};
					};
				};
				#print "new:\n $new_coordspecial\n";

			if($new_coordspecial eq ''){
					
				print OUTC "$new_coord[0]";
				print OUTC "$new_coord[1]";
				print OUTC "$new_coord[2]";
				$new_coord[3]=sprintf"%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$longa,$longb2;
				print OUTC "$new_coord[3]";
				print OUTC "$new_coord[4]";
				foreach $rf (0..@addh-1){
					print OUTC $addh[$rf];
					$differentielx++ if($addh[$rf]=~/ X /);
				};
				print OUTC "$new_coord[5]";
				
				$differentiel=$differentiel+$longb2;
				foreach $rf (0..@addh-1){
					print OUTC $addh_bond[$rf];
				};
				
				if($rgnb != 0){
	                        	$rgx=sprintf"%5.3f",$rgx/$rgnb;
        	                	$rgy=sprintf"%5.3f",$rgy/$rgnb;
                	        	$rgz=sprintf"%5.3f",$rgz/$rgnb;
				};

				print OUTC "M  END\n";
				print OUTC "> <ncycles>\n";
				print OUTC "1\n";
				print OUTC "\n";
				print OUTC "> <natomescycle>\n";
				print OUTC "$long\n";
				print OUTC "\n";
				print OUTC "> <cyclearom>\n";
				print OUTC "$nar\n";
				print OUTC "\n";
				print OUTC "> <nhetat>\n";
				print OUTC "$nhetat\n";
				print OUTC "\n";
	
				if($nbatomar == $long){
        	                	print OUTC "> <ARcenter>\n";
	                        	print OUTC "$rgx $rgy $rgz\n";
                	        	print OUTC "\n";
				}
				else{
					print OUTC "> <LIPcenter>\n";
					print OUTC "$rgx $rgy $rgz\n";
					print OUTC "\n";
				};

				if($point ne ''){
					print OUTC "> <POINTS>\n";
					print OUTC "$point\n";
					print OUTC "\n";
					print OUTC "> <ORIGIN>\n";
					print OUTC "$pointo\n";
					print OUTC "\n";
				};
				print OUTC "$lignedata" if($lignedata ne "");
				print OUTC "\$\$\$\$\n";
			}
			else{
				print OUTD "$new_coord[0]";
				print OUTD "$new_coord[1]";
				print OUTD "$new_coord[2]";
				$new_coord[3]=sprintf"%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$longa,$longb2;
				print OUTD "$new_coord[3]";
				print OUTD "$new_coord[4]";
				foreach $rf (0..@addh-1){
				 	print OUTD $addh[$rf];
					$differentielx++ if($addh[$rf]=~/ X /);
				};
				print OUTD "$new_coord[5]";
				print OUTD "$new_coordspecial" if($new_coordspecial ne "");
				
		       		$differentiel=$differentiel+$longb2;
				foreach $rf (0..@addh-1){
					print OUTD $addh_bond[$rf];
				};

				$rgx=sprintf"%5.3f",$rgx/$rgnb;
				$rgy=sprintf"%5.3f",$rgy/$rgnb;
				$rgz=sprintf"%5.3f",$rgz/$rgnb;

				print OUTD "M  END\n";
				print OUTD "> <ncycles>\n";
				print OUTD "$longgetstar\n";
				print OUTD "\n";
				print OUTD "> <natomescycle>\n";
				print OUTD "$specialib\n";
				print OUTD "\n";
				print OUTD "> <cyclearom>\n";
				print OUTD "$nar\n";
				print OUTD "\n";
				print OUTD "> <nhetat>\n";
				print OUTD "$nhetat\n";
				print OUTD "\n";
				if($nbringar > 0){
					print OUTD "> <nbringar>\n";
					print OUTD "$nbringar\n";
					print OUTD "\n";
				};	
			
				if($nbatomar == $specialib){
					print OUTD "> <ARcenter>\n";
					print OUTD "$rgx $rgy $rgz\n";
					print OUTD "\n";
				}
				else{
					print OUTD "> <LIPcenter>\n";
					print OUTD "$rgx $rgy $rgz\n";
					print OUTD "\n";
				};

				if($point ne ''){
					print OUTD "> <POINTS>\n";
					print OUTD "$point\n";
					print OUTD "\n";
					print OUTD "> <ORIGIN>\n";
					print OUTD "$pointo\n";
					print OUTD "\n";
				};
				print OUTD "$lignedata" if($lignedata ne "");
				print OUTD "\$\$\$\$\n";
				
			};	
				
			};
		};

close(OUTC);



###################################################################### Linker (>= 3 atomes) en SDF
#print "@copylign\n";

open(OUTC,">>linker_tmp.sdf");
	
	$cyclepoint=" ";
	foreach $k (1..@copylign-1){
		if ($copylign[$k] eq ''){
			$cyclepoint=$cyclepoint." $k ";
		};
	};


foreach $k (1..@copylign-1){
        $nhetat=0;
	$point="";
	@corespond='';
	@addh='';
	@addh_bond='';
	$addhi=0;
	$longa=0;
	$longb=0;
	@atom="";
	$nrot=0;
	$suitehalo=" H F Cl Br I ";
	$suiteheavy=" C N O S ";
	$rgx=0;
	$rgy=0;
	$rgz=0;
	$rgnb=0;
	if ($copylign[$k] ne ''){
		$pile=" $k ";
		@get8=split(' ',$pile);
		
		@get7 = split(' ',$lignecor[$get8[0]]);
		
		$nhetat++ if (($get7[3]=~/O/)||($get7[3]=~/S/)||($get7[3]=~/N/));
		$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get7[0],$get7[1],$get7[2],$blanc,$get7[3],$blanc,$get7[4],$get7[5],$get7[6];	
		$atom[$addhi+1]=$get7[3];

		# Geometric center
		if($lignecor[$get8[0]]!~/ H / && $lignecor[$get8[0]]!~/ X /){
			$rgx=$rgx+$get7[0];
			$rgy=$rgy+$get7[1];
			$rgz=$rgz+$get7[2];
			$rgnb++;
		};

		$corespond[$get8[0]]=$addhi+1;	
		$addhi++;
		$longa++;
		$copylign[$get8[0]]='';
			
		while($get8[0] ne ''){
			
			@get = split(' ',$listb[$get8[0]]);
			@get2 = split(' ',$typeb[$get8[0]]);

			foreach $l (0..@get-1){
				if($copylign[$get[$l]] ne ''){
					
					$pile=$pile." ".$get[$l]." ";
					@get7 = split(' ',$lignecor[$get[$l]]);
					$nhetat++ if (($get7[3]=~/O/)||($get7[3]=~/S/)||($get7[3]=~/N/));
		        		$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get7[0],$get7[1],$get7[2],$blanc,$get7[3],$blanc,$get7[4],$get7[5],$get7[6];	
					$atom[$addhi+1]=$get7[3];
					$y=$longa+1;
					$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$corespond[$get8[0]],$y,$get2[$l];	
		
					# Geometric center
					if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
						$rgx=$rgx+$get7[0];
						$rgy=$rgy+$get7[1];
						$rgz=$rgz+$get7[2];
						$rgnb++;
					};
#print "$suitehalo <=> $atom[$y]\n";

					#rotatable (single bond, not in ring, not terminal)
					if($get2[$l] == 1 && $atom[$corespond[$get8[0]]] ne "" && $atom[$y] ne "" && $suitehalo!~/ $atom[$corespond[$get8[0]]] / && $suitehalo!~/ $atom[$y] /){
						
						# search if $get8[0] non terminal heavy atom
						$nrottmp1=0;
						@narp=split(' ',$listb[$get8[0]]);
						foreach $narpi (0..@narp-1){
							if($narp[$narpi] != $get[$l]){
								#print "\t$narp[$narpi]\n";
								if($suiteheavy=~/ $type[$narp[$narpi]] /){
									$nrottmp1=1;
									#print "ok1\n";
								};
							};
						};	
					       	
						# search if $get[$l] non terminal heavy atom
						@narp=split(' ',$listb[$get[$l]]);
						$nrottmp=0;
						foreach $narpi (0..@narp-1){
							if($narp[$narpi] != $get8[0]){
								#print "\t$narp[$narpi]\n";
								if($suiteheavy=~/ $type[$narp[$narpi]] /){
									$nrottmp=1;
									#print "ok\n";
								};
							};
						};	
						
						$nrot++ if($nrottmp && $nrottmp1);
						#print "in $get[$l]-$get8[0] $get[$l] heavy not terminal\n" if($nrottmp && $nrottmp1);	
					}; 								

					$corespond[$get[$l]]=$addhi+1;	
					$addhi++;
					$longa++;
					$longb++;
					$copylign[$get[$l]]='';
				}
				elsif($copylign[$get[$l]] eq '' && $cyclepoint=~/ $get[$l] /){
					@get7 = split(' ',$lignecor[$get[$l]]);
		        		$get7[3]="X";
		        		$nhetat++ if (($get7[3]=~/O/)||($get7[3]=~/S/)||($get7[3]=~/N/));
		        		$addh[$addhi]=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get7[0],$get7[1],$get7[2],$blanc,$get7[3],$blanc,$get7[4],$get7[5],$get7[6];	
					$atom[$addhi+1]=$get7[3];
					$y=$longa+1;
					$addh_bond[$addhi]=sprintf"%3s%3s%3s  0  0  0  0\n",$corespond[$get8[0]],$y,1;	


					# Geometric center
					if($addh[$addhi]!~/ H / && $addh[$addhi]!~/ X /){
						$rgx=$rgx+$get7[0];
						$rgy=$rgy+$get7[1];
						$rgz=$rgz+$get7[2];
						$rgnb++;
					};

					#rotatable
					if($atom[$corespond[$get8[0]]] ne "" && $atom[$y] ne "" && $suitehalo!~/ $atom[$corespond[$get8[0]]] / && $suitehalo!~/ $atom[$y] /){
						# search if $get8[0] non terminal heavy atom
						$nrottmp1=0;
						@narp=split(' ',$listb[$get8[0]]);
						foreach $narpi (0..@narp-1){
							if($narp[$narpi] != $get[$l]){
								if($suiteheavy=~/ $type[$narp[$narpi]] /){
									$nrottmp1=1;
								};
							};
						};	
					       	
						# search if $get[$l] non terminal heavy atom
						# It is since $get[$l] belongs to a ring system
						$nrottmp=1;
						
						$nrot++ if($nrottmp && $nrottmp1);
						#print "in $get[$l]-$get8[0] $get[$l] in ring heavy not terminal\n" if($nrottmp && $nrottmp1);
						
					}
					elsif($atom[$corespond[$get8[0]]] eq "" || $atom[$y] eq ""){
						print "nrotatable: atom no $corespond[$get8[0]] or atom no $atom[$y] unknown !\n"; 
					};

					$addhi++;
					$longa++;
					$longb++;
					if($point eq ''){
						$point=$corespond[$get8[0]];
						$pointo=$get8[0]."_".$get[$l];
					}
					else{
						$point=$point."-".$corespond[$get8[0]];
						$pointo=$pointo."-".$get8[0]."_".$get[$l];
					};
				};
			};
			$pile=~s/ $get8[0] //g;
			@get8=split(' ',$pile);
			#print "pile = $pile\n";
		};
	};
	
	if($addhi > 0){
		print OUTC "\n";
		print OUTC "  -ISIS-  08029913582D\n";
		print OUTC "\n";
		printf OUTC "%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$longa,$longb;
		foreach $rf (0..@addh-1){
			print OUTC $addh[$rf];
			$differentielx++ if($addh[$rf]=~/ X /);
		};
				
		foreach $rf (0..@addh-1){
			print OUTC $addh_bond[$rf];
		};
		$differentiel=$differentiel+$longb;

		if($rgnb != 0){ 	
		$rgx=sprintf"%5.3f",$rgx/$rgnb;
		$rgy=sprintf"%5.3f",$rgy/$rgnb;		
		$rgz=sprintf"%5.3f",$rgz/$rgnb;
		};
		
		print OUTC "M  END\n";
		print OUTC "> <nhetat>\n";
		print OUTC "$nhetat\n";
		print OUTC "\n";
		print OUTC "> <natomes>\n";
		print OUTC "$longa\n";
		print OUTC "\n";
		print OUTC "> <nrotatable>\n";
		print OUTC "$nrot\n";
		print OUTC "\n";	
		print OUTC "> <LIPcenter>\n";
		print OUTC "$rgx $rgy $rgz\n";
		print OUTC "\n";
		if($point ne ''){
			print OUTC "> <POINTS>\n";
			print OUTC "$point\n";
			print OUTC "\n";
			print OUTC "> <ORIGIN>\n";
			print OUTC "$pointo\n";
			print OUTC "\n";
		};
		print OUTC "$lignedata" if($lignedata ne "");
		print OUTC "\$\$\$\$\n";	
		
	};
};
	

	
close(OUTC);
	
######################################################################

};

close(OUTC2);


########## ########## ########## ########## ########## ########## ##########
########## ########## ########## ########## ########## ########## ##########
########## ########## ########## ########## ########## ########## ##########



				 				
########## ########## ########## ########## ########## ########## ##########
########## ########## ########## ########## ########## ########## ##########
########## ########## ########## ########## ########## ########## ##########

sub fusion{
	local($voisin)=@_;


			@get2= split(' ',$cyclef[$voisin]);
			$long2=@get2;

#print "avant @get2\n";
			@get3= split(' ',$cyclef[$voisin]);
			$long3=@get3;

			@get3n='';
#print "put @put\n";
			#print "der $ligneb[$lb-1]\n";
			
			foreach $rg (0..@put-1){
				#print "$put[$rg]\n";
				foreach $rg2 (0..@get2-1){
					if ($put[$rg] == $get2[$rg2]){
						#print "$put[$rg] = $get2[$rg2]\n";
						# verifie que le bond est fait entre le i et le i+1
						foreach $rg3 (0..@put-1){
							#print "$put[$rg3] / $get2[$rg2+1]\n";
							
							# regarde liaison d'avant
							if($get2[$rg2] == $get2[@get2-1]){
								if ($put[$rg3] == $get2[0]){
									#print "elements $get2[$rg2] et $get2[0] deja vus \n";

									@getb = split(' ',$listb[$get2[$rg2]]);
									foreach $rfb (0..@getb-1){
										$checkbi=$rfb if ($get2[0] == $getb[$rfb]);
									};
									@getb = split(' ',$typeb[$get2[$rg2]]);

									$checkb1=sprintf"%3s%3s%3s  0  0  0  0",$rg3+1,$rg+1,$getb[$checkbi];
									$checkb2=sprintf"%3s%3s%3s  0  0  0  0",$rg+1,$rg3+1,$getb[$checkbi];

									$dejavvucoupleb=0;
									foreach $rfb (0..@ligneb-1){
										if($ligneb[$rfb] eq $checkb1 || $ligneb[$rfb] eq $checkb2){
											#print "couple-bond deja vu\n";
											$dejavvucoupleb=1;
										};
									};
									if($dejavvucoupleb==0){
										#print "couple-bond a ajouter !\n";
										$ligneb[$lb]=$checkb1;
										$lb++;
									};

								};
							}
							else{
								if ($put[$rg3] == $get2[$rg2+1]){
									#print "elements $get2[$rg2] et $get2[$rg2+1] deja vus \n";
									
									@getb = split(' ',$listb[$get2[$rg2]]);
									foreach $rfb (0..@getb-1){
										$checkbi=$rfb if ($get2[$rg2+1] == $getb[$rfb]);
									};
									@getb = split(' ',$typeb[$get2[$rg2]]);

									$checkb1=sprintf"%3s%3s%3s  0  0  0  0",$rg3+1,$rg+1,$getb[$checkbi];
									$checkb2=sprintf"%3s%3s%3s  0  0  0  0",$rg+1,$rg3+1,$getb[$checkbi];
#print "$checkb1\n";


									$dejavvucoupleb=0;
									foreach $rfb (0..@ligneb-1){
										if($ligneb[$rfb] eq $checkb1 || $ligneb[$rfb] eq $checkb2){
											#print "couple-bond deja vu\n";
											$dejavvucoupleb=1;
										};
									};
									if($dejavvucoupleb==0){
										#print "couple-bond a ajouter !\n";
										$ligneb[$lb]=$checkb1;
										$lb++;
									};

								};
							};
							
							# regarde liaison d'apres
							if($get2[$rg2] == $get2[0]){
								if ($put[$rg3] == $get2[@get2-1]){
									#print "elements $get2[$rg2] et $get2[@get2-1] deja vus \n";

									@getb = split(' ',$listb[$get2[$rg2]]);
									foreach $rfb (0..@getb-1){
										$checkbi=$rfb if ($get2[@get2-1] == $getb[$rfb]);
									};
									@getb = split(' ',$typeb[$get2[$rg2]]);

									$checkb1=sprintf"%3s%3s%3s  0  0  0  0",$rg3+1,$rg+1,$getb[$checkbi];
									$checkb2=sprintf"%3s%3s%3s  0  0  0  0",$rg+1,$rg3+1,$getb[$checkbi];

									$dejavvucoupleb=0;
									foreach $rfb (0..@ligneb-1){
										if($ligneb[$rfb] eq $checkb1 || $ligneb[$rfb] eq $checkb2){
											#print "couple-bond deja vu\n";
											$dejavvucoupleb=1;
										};
									};
									if($dejavvucoupleb==0){
										#print "couple-bond a ajouter !\n";
										$ligneb[$lb]=$checkb1;
										$lb++;
									};

								};
							}
							else{
								if ($put[$rg3] == $get2[$rg2-1]){
									#print "elements $get2[$rg2] et $get2[$rg2+1] deja vus \n";
									
									@getb = split(' ',$listb[$get2[$rg2]]);
									foreach $rfb (0..@getb-1){
										$checkbi=$rfb if ($get2[$rg2-1] == $getb[$rfb]);
									};
									@getb = split(' ',$typeb[$get2[$rg2]]);

									$checkb1=sprintf"%3s%3s%3s  0  0  0  0",$rg3+1,$rg+1,$getb[$checkbi];
									$checkb2=sprintf"%3s%3s%3s  0  0  0  0",$rg+1,$rg3+1,$getb[$checkbi];
#print "$checkb1\n";


									$dejavvucoupleb=0;
									foreach $rfb (0..@ligneb-1){
										if($ligneb[$rfb] eq $checkb1 || $ligneb[$rfb] eq $checkb2){
											#print "couple-bond deja vu\n";
											$dejavvucoupleb=1;
										};
									};
									if($dejavvucoupleb==0){
										#print "couple-bond a ajouter !\n";
										$ligneb[$lb]=$checkb1;
										$lb++;
									};

								};
							};
						};

						$get2[$rg2]='';
					};
				};
			};
#print "apres elim @get2\n";

			$rh=@put;
			foreach $rg (0..@get2-1){

				if ($get2[$rg] ne ''){
					$put[$rh]=$get2[$rg];

					$narp=" ".$typeb[$put[$rh]]." ";
					$nocc1=$narp=~s/1/1/g;
					$nocc2=$narp=~s/2/2/g;
					if($type[$put[$rh]] eq 'C'){
						$nbatomar++ if($nocc1==2 && $nocc2==1);
					}
					elsif($type[$put[$rh]] eq 'O'){
						$nbatomar++ if($nocc1==2 && $long2==5);
					}
					elsif($type[$put[$rh]] eq 'N'){
						if($nocc1==1 && $nocc2==1){
							$nbatomar++;
						}
						elsif($nocc1==3 && $long2==5){
							$nbatomar++;
						};
					}
					elsif($type[$put[$rh]] eq 'S'){
						$nbatomar++ if($nocc1==2 && $long2==5);
					};
					
					$get3n[$rg]=$rh+1;
					$rh++;
				}
				else{
					foreach $rg2 (0..@put-1){
						if ($get3[$rg] eq $put[$rg2]){
							$get3n[$rg]=$rg2+1;
						};
					};
				};
			};


			foreach $rg (0..@get3-1){
				if (($get2[$rg] eq '') && ($get2[$rg+1] eq '')){
					$situ=10000;
				}
				else{
					$situ=10000;
					@getb = split(' ',$listb[$get3[$rg]]);
					foreach $rfb (0..@getb-1){
						$situ=$rfb if ($get3[$rg+1] == $getb[$rfb]);
						if ($get2[$rg] ne ''){
							@getb2 = split(' ',$listb[$getb[$rfb]]);
							if (@getb2 == 1){
								$put2[$p2]=$getb[$rfb];
								@get3b = split(' ',$typeb[$get3[$rg]]);
								$x=$get3n[$rg];
								$y=$p2;
								$z="  0  0  0  0";
								$ligneb2[$lb2]=sprintf"%3s%3s%3s  0  0  0  0",$x,$y,$get3b[$rfb];
								$lb2++;
								$p2++;
							}
							else{
								$flagb=1;
								foreach $b2 (0..@getb2-1){
									@getb3 = split(' ',$listb[$getb2[$b2]]);
									if ($getb2[$b2] != $get3[$rg]){
										$flagb=0 if (@getb3 > 1);
									};
								};
								if ($flagb){
									$put2[$p2]=$getb[$rfb];
									@get3b = split(' ',$typeb[$get3[$rg]]);
									$x=$get3n[$rg];
									$y=$p2;
									$z="  0  0  0  0";
									$ligneb2[$lb2]=sprintf"%3s%3s%3s  0  0  0  0",$x,$y,$get3b[$rfb];
									$lb2++;
									$p2++;
								
									foreach $b2 (0..@getb2-1){
										if ($getb2[$b2] != $get3[$rg]){
											$put3[$p3]=$getb2[$b2];
											@get4b = split(' ',$typeb[$getb[$rfb]]);
											$x=$p2-1;
											$y=$p3;
											$z="  0  0  0  0";
											$ligneb3[$lb3]=sprintf"%3s%3s%3s  0  0  0  0",$x,$y,$get4b[$b2];
											$lb3++;
											$p3++;
										};
									};
								};
							};
						};
					};
					if ($situ != 10000){
						@get3b = split(' ',$typeb[$get3[$rg]]);
						$x=$get3n[$rg];
						$y=$get3n[$rg+1];
						$z="  0  0  0  0";
						$ligneb[$lb]=sprintf"%3s%3s%3s  0  0  0  0",$x,$y,$get3b[$situ];
						#print "ici3 $ligneb[$lb]\n";
						$lb++;
					};
				};
			};
			if (($get2[0] eq '') && ($get2[@get2-1] eq '')){
				$situ=10000;
			}
			else{
				$situ=10000;
				@getb = split(' ',$listb[$get3[0]]);
				foreach $rfb (0..@getb-1){
					$situ=$rfb if ($get3[@get3-1] == $getb[$rfb]);
				};
				if ($situ != 10000){
					@get3b = split(' ',$typeb[$get3[0]]);
					$x=$get3n[0];
					$y=$get3n[@get2-1];
					$z="  0  0  0  0";
					$ligneb[$lb]=sprintf"%3s%3s%3s  0  0  0  0",$x,$y,$get3b[$situ];
					$lb++;
				};
			};

};

