#!/usr/bin/perl

my($file)=@ARGV;

if($file eq '' ){
print "usage: \n nosel <file.sdf>\n";
exit;
};
		
	unlink "one.sdf" if(-e "one.sdf");
	unlink "two.sdf" if(-e "two.sdf");
	unlink "dissociated.sdf" if(-e "dissociated.sdf");
	open(IN,">one.sdf");
	close(IN);
	open(IN,">dissociated.sdf");
	close(IN);
	
	$flagnew=1;
	$moli=0;
	open(IN,"<$file");
	while(<IN>){
		if($flagnew){
			$masse=0;
			$moli++;
			$compt=0;
			$ig=1;
			$jg=0;
			@strx='';
			@stry='';
			@strz='';
			@atom='';
			@coval='';
			@fonc='';
			@ifonc='';
			@covfonc='';
			@bond='';
			@listb='';
			@typeb='';
			$blanc=' ';	
			@radius='';
			@lignebond='';		
			$atomlourd=0;
			$flagnew=0;
			
			@sdf='';
		};
		
		@getstr = split(' ',$_);
		
		$sdf[$compt]=$_;
		$compt++;

		if (($compt > 4) && ($ig <= $istratom)){
			$strx[$ig]=$getstr[0];
			$stry[$ig]=$getstr[1];
			$strz[$ig]=$getstr[2];
			$atom[$ig]=$getstr[3];
			$atomlourd ++ if($getstr[3] ne 'H' && $getstr[3] ne 'X');
			$radius[$ig]=$tabR{$getstr[3]};
			$masse=$masse+$tabmm{$getstr[3]};
			$ig++;
		};
		if (($compt > 4) && ($ig > $istratom) && ($jg <=$istrbond)){
			if ($jg == 0){
				$jg++;
			}
			else{

                                @coller=split(' *',$getstr[0]);
                                @coller2=split(' *',$getstr[1]);
                                if(@coller==6 && $getstr[1] ne ""){
                                        $getstr[0]=$coller[0].$coller[1].$coller[2];
                                        $getstr[2]=$getstr[1];
                                        $getstr[1]=$coller[3].$coller[4].$coller[5];
                                }
                                elsif(@coller==6 && $getstr[1] eq ""){
                                        $getstr[0]=$coller[0].$coller[1];
                                        $getstr[1]=$coller[2].$coller[3].$coller[4];
                                        $getstr[2]=$coller[5];
                                }
                                elsif(@coller==5){
                                        if($_=~/^\s/){
                                                $getstr[0]=$coller[0].$coller[1];
                                                $getstr[2]=$getstr[1];
                                                $getstr[1]=$coller[2].$coller[3].$coller[4];
                                        }
                                        else{
                                                $getstr[0]=$coller[0].$coller[1].$coller[2];
                                                $getstr[2]=$getstr[1];
                                                $getstr[1]=$coller[3].$coller[4];
                                        };
                                }
                                elsif(@coller==4){
                                        if($_=~/^\s/){
                                                $getstr[0]=$coller[0];
                                                $getstr[2]=$getstr[1];
                                                $getstr[1]=$coller[1].$coller[2].$coller[3];
                                        }
                                        else{
                                                $getstr[0]=$coller[0].$coller[1].$coller[2];
						$getstr[2]=$getstr[1];
                                                $getstr[1]=$coller[3];
                                        };					
                                }
                                elsif(@coller2==4){
                                        $getstr[1]=$coller2[0].$coller2[1].$coller2[2];
                                        $getstr[2]=$coller2[3];
                                }
                                elsif(@coller==7){
                                        $getstr[0]=$coller[0].$coller[1].$coller[2];
                                        $getstr[1]=$coller[3].$coller[4].$coller[5];
                                        $getstr[2]=$coller[6];
                                };

				$bond[$getstr[0]]=$bond[$getstr[0]].$blanc.$getstr[1].$blanc.$getstr[2];
				$listb[$getstr[0]]=$listb[$getstr[0]].$blanc.$getstr[1];
				$typeb[$getstr[0]]=$typeb[$getstr[0]].$blanc.$getstr[2];

				$bond[$getstr[1]]=$bond[$getstr[1]].$blanc.$getstr[0].$blanc.$getstr[2];
				$listb[$getstr[1]]=$listb[$getstr[1]].$blanc.$getstr[0];
				$typeb[$getstr[1]]=$typeb[$getstr[1]].$blanc.$getstr[2];


				$fonc[$getstr[0]]=$fonc[$getstr[0]].$blanc.$getstr[2].'-'.$atom[$getstr[1]].$blanc;
				$ifonc[$getstr[0]]=$ifonc[$getstr[0]].$blanc.$getstr[1].$blanc;
				$covfonc[$getstr[0]]=$covfonc[$getstr[0]].$blanc.$getstr[2];
				$coval[$getstr[0]]=$coval[$getstr[0]]+$getstr[2];

				$fonc[$getstr[1]]=$fonc[$getstr[1]].$blanc.$getstr[2].'-'.$atom[$getstr[0]].$blanc;
				$ifonc[$getstr[1]]=$ifonc[$getstr[1]].$blanc.$getstr[0].$blanc;
				$covfonc[$getstr[1]]=$covfonc[$getstr[1]].$blanc.$getstr[2];
				$coval[$getstr[1]]=$coval[$getstr[1]]+$getstr[2];
				$lignebond[$jg]=$_;
				$jg++;
			};
		};
		if ($compt == 4){
			$istratom=$getstr[0];
			$istrbond=$getstr[1];

			@coller=split(' *',$istratom);
			if(@coller>3 && @coller==6){
				$istratom=$coller[0].$coller[1].$coller[2];
				$istrbond=$coller[3].$coller[4].$coller[5];
			}
			elsif(@coller>3 && @coller==5){
				if($_=~/^\s/){
					$istratom=$coller[0].$coller[1];
					$istrbond=$coller[2].$coller[3].$coller[4];
				}
				else{
					$istratom=$coller[0].$coller[1].$coller[2];
					$istrbond=$coller[3].$coller[4];
				};
			};

			#$compt++;
		};
		if ($_=~/\$\$\$\$/){
			$flagnew=1;
			if($istratom > 0){
			&convert;
			
			if($one){
			        print "Mol$moli 1 ($istratom atomes)\n";
				open(OUT,">>one.sdf");
				foreach $p (0..@sdf-1){
					print OUT "$sdf[$p]";
				};
				close(OUT);
			}
			else{
				print "Mol$moli 2 --> two.sdf ($istratom atomes)\n";
				open(OUT,">>two.sdf");
				foreach $p (0..@sdf-1){
					print OUT "$sdf[$p]";
				};
				close(OUT);
			};
			};
		};
		
		
	};
	close(IN);
        print "dissociated.sdf (>= 5 heavy atoms) contains separated parts of salt - one.sdf contains sdf that do not have salt - two.sdf contains sdf of salts\n";
	print "Recommended: 'cat one.sdf dissociated.sdf > new_file.sdf'\n";

###################################################################################
###################################################################################


sub convert{
	
	$one=0;
	$suite=" 1 ";
	$continue=1;
	$bi=1;
	$pos=0;
	while($continue){
                @get=split(' ',$ifonc[$bi]);
                foreach $k (0..@get-1){	
                	$suite=$suite."$get[$k] " if ($suite!~/ $get[$k] /); 	
                };	
                @get2=split(' ',$suite);
                $pos++;
                if($pos == @get2){
                	 $continue=0;
                }
                else{
                	$bi=$get2[$pos];
                };
	};
	
	#print "$suite\n";
        @get2=split(' ',$suite);
        $nbacc= @get2;
        $acc=$suite;

        if(@get2 == $istratom){
        	$one=1;
        }
        else{
        	#SPLIT LES FICHIERS
 		  	
		@geto=split(' ',$suite);
		$atl=0;
		foreach $m (0..@geto-1){
			$atl++ if($atom[$geto[$m]] ne 'H' && $atom[$geto[$m]] ne 'X');
		};       	
		&print($suite) if($atl > 4);
		
		$continue=1;
		while($continue){
		        $suite=" ";
		        $bi='';
		        foreach $k (1..$istratom){	
		        	$bi=$k if($acc!~/ $k /);
		        };
	                $suite=$suite."$bi ";
	
	                $cont=1;
	                $pos=0;
	                while($cont){
		        	@get=split(' ',$ifonc[$bi]);
		        	foreach $k (0..@get-1){	
		               		$suite=$suite."$get[$k] " if ($suite!~/ $get[$k] /); 	
		        	};
		        	@get2=split(' ',$suite);
		        	$pos++;
                		if($pos == @get2){
                	 		$cont=0;
                		}
               	 		else{
                			$bi=$get2[$pos];
                		};		        	
		        };
		        
			$acc=$acc.' '.$suite;
			@getd=split(' ',$acc);
		  	$continue=0 if(@getd >= $istratom);
		  	$continue=0 if($bi eq '');
		  	
		  	@geto=split(' ',$suite);
		  	$atl=0;
		  	foreach $m (0..@geto-1){
		  		$atl++ if($atom[$geto[$m]] ne 'H' && $atom[$geto[$m]] ne 'X');
		  	};
			&print($suite) if($atl > 4);
		};
        };



        $one;
};

######################################################################
######################################################################

sub print{
         local($atomes)=@_;
#print "$atomes\n";

        $atomes=' '.$atomes.' ';
	@getp=split(' ',$atomes);
	$longa=@getp;
	
	$longb=0;
	$gh=1;
	@ligneb='';
	foreach $f (($istratom+5-1)..@sdf-1){ # 5-1 to begin index 0 see $compt++
		#print "$f $sdf[$f]\n";
	        @getb=split(' ',$sdf[$f]);
	
                               @coller=split(' *',$getb[0]);
                                @coller2=split(' *',$getb[1]);
                                if(@coller==6 && $getb[1] ne ""){
                                        $getb[0]=$coller[0].$coller[1].$coller[2];
                                        $getb[2]=$getb[1];
                                        $getb[1]=$coller[3].$coller[4].$coller[5];
                                }
                                elsif(@coller==6 && $getb[1] eq ""){
                                        $getb[0]=$coller[0].$coller[1];
                                        $getb[1]=$coller[2].$coller[3].$coller[4];
                                        $getb[2]=$coller[5];
                                }
                                elsif(@coller==5){
                                        if($sdf[$f]=~/^\s/){
                                                $getb[0]=$coller[0].$coller[1];
                                                $getb[2]=$getb[1];
                                                $getb[1]=$coller[2].$coller[3].$coller[4];
                                        }
                                        else{
                                                $getb[0]=$coller[0].$coller[1].$coller[2];
                                                $getb[2]=$getb[1];
                                                $getb[1]=$coller[3].$coller[4];
                                        };
                                }
                                elsif(@coller==4){
                                        if($sdf[$f]=~/^\s/){
                                                $getb[0]=$coller[0];
                                                $getb[2]=$getb[1];
                                                $getb[1]=$coller[1].$coller[2].$coller[3];
                                        }
                                        else{
                                                $getb[0]=$coller[0].$coller[1].$coller[2];
						$getb[2]=$getb[1];
                                                $getb[1]=$coller[3];
                                        };					
                                }
                                elsif(@coller2==4){
                                        $getb[1]=$coller2[0].$coller2[1].$coller2[2];
                                        $getb[2]=$coller2[3];
                                }
                                elsif(@coller==7){
                                        $getb[0]=$coller[0].$coller[1].$coller[2];
                                        $getb[1]=$coller[3].$coller[4].$coller[5];
                                        $getb[2]=$coller[6];
                                };

	
		$gh=0 if($sdf[$f]=~/^>/ || $getb[0] eq '' || $sdf[$f]=~/^M/ || $sdf[$f]=~/^\$\$\$\$/);
	        if($gh){
			#print "getb0 $getb[0]\n";
	        	if($atomes=~/ $getb[0] /){
				#print "\t $getb[0]\n";
	        		$a1='';
	        		$a2='';
        			foreach $p (0..@getp-1){
        				$a1=($p+1) if($getb[0] eq $getp[$p]);
        				$a2=($p+1) if($getb[1] eq $getp[$p]);
        	                };
        	                $z="  0  0  0  0";
				#$ligneb[$longb]=pack("A3A3A3A12",$a1,$a2,$getb[2],$z);
        	                $ligneb[$longb]=sprintf"%3s%3s%3s  0  0  0  0",$a1,$a2,$getb[2];
				#	print "$ligneb[$longb]\n";
				$longb++;
        		};
        	};
 	};


	open(OUT,">>dissociated.sdf");
	foreach $p (0..2){
		print OUT "$sdf[$p]";
	};
	
	printf OUT "%3s%3s  0  0  0  0  0  0  0  0999 V2000\n",$longa,$longb;
	
 	foreach $f (0..@getp-1){
		#print "$getp[$f] \t $sdf[$getp[$f]+3]";
               print OUT "$sdf[$getp[$f]+3]";
  	};

	foreach $f (0..@ligneb-1){
               print OUT "$ligneb[$f]\n";
  	};	
			
	print OUT "M  END\n";
	
	$ecritfin=0;
	foreach $mp (0..@sdf-1){
		$ecritfin=1 if($sdf[$mp]=~/^>/);	
		print OUT "$sdf[$mp]" if($ecritfin);
	};
	print OUT "\$\$\$\$\n" if($ecritfin==0);	
	close(OUT);		
};

