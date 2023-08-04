#!/usr/bin/perl

	my($filesdf,$fsdfo)=@ARGV;
	if($filesdf eq '' || $fsdfo eq ''){
	        print "usage: \n .pl <file.sdf> <name of output file>\nreturn a file with 3D molecules only (having non null z coordinates)\n";
	        exit;
	};
	
	$moli=0;
	$flagnew=1;
	open(MPO,">$fsdfo");
	open(MOL,"<$filesdf");
	while(<MOL>){
		if($flagnew){
			$masse=0;
			$nbx=0;
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
			$blanc=' ';	
			@radius='';
			$atomlourd=0;
			$flagnew=0;
			@line='';
		};
		@getstr = split(' ',$_);
		$line[$compt]=$_;
		$compt++;
		if (($compt > 4) && ($ig <= $istratom)){
			$strx[$ig]=$getstr[0];
			$stry[$ig]=$getstr[1];
			$strz[$ig]=$getstr[2];
			$atom[$ig]=$getstr[3];
			$atomlourd ++ if($getstr[3] ne 'H');
			$radius[$ig]=$tabR{$getstr[3]};
			$masse=$masse+$tabmm{$getstr[3]};
			#$nbx++ if($getstr[3] eq 'X');
			$warn2d++ if($strz[$ig]==0);
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

				$fonc[$getstr[0]]=$fonc[$getstr[0]].$blanc.$getstr[2].'-'.$atom[$getstr[1]].$blanc;
				$ifonc[$getstr[0]]=$ifonc[$getstr[0]].$blanc.$getstr[1].$blanc;
				$covfonc[$getstr[0]]=$covfonc[$getstr[0]].$blanc.$getstr[2];
				$coval[$getstr[0]]=$coval[$getstr[0]]+$getstr[2];

				$fonc[$getstr[1]]=$fonc[$getstr[1]].$blanc.$getstr[2].'-'.$atom[$getstr[0]].$blanc;
				$ifonc[$getstr[1]]=$ifonc[$getstr[1]].$blanc.$getstr[0].$blanc;
				$covfonc[$getstr[1]]=$covfonc[$getstr[1]].$blanc.$getstr[2];
				$coval[$getstr[1]]=$coval[$getstr[1]]+$getstr[2];

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
			$moli++;
			$flagnew=1;
			$masse=sprintf"%4.2f",$masse;
			#print "$masse\n";
			#print "$nbx\n";
			if($warn2d < $istratom){
				#print "$moli in 3D\n";
				foreach $k (0..@line-1){
					print MPO $line[$k];
				};
			};
		};
		
	};
	close(MOL);
	close(MPO);
	
#################################################################################################"
#################################################################################################"


