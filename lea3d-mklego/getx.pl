#!/usr/bin/perl


$f='';
my($f)=@ARGV;

if($f eq ''){
	print "usage: getx <file.sdf>\n";
	exit;
};
$ecrit=1;
$point="";
$point2="";
$pasvu=1;
$comp=0;
$file="new_".$f;	
$i=1;
$j=0;
open(DOC,">$file");
open(IN,"<$f");
while(<IN>){

        @get = split(' ',$_);
        $comp++;

        if ($comp == 4){
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
			
	};
	if (($comp > 4) && ($i <= $nbatom)){
		$comp2=$comp-4;
		
        	$point=$point."-".$comp2 if ($get[3] eq "X");
        	
        	#print "point $point $get[3]\n" if ($get[3] eq "X");
        	
        	$i++;
	};
	
	if (($comp > 4) && ($i > $nbatom) && ($j <= $nbond)){
	 	if ($j == 0){
			$j++;
			$point=$point."-";
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
	
				#print "point $point & 1= $get[1] \n";
				#print "point $point & 2= $get[2] \n";
				
				if($point=~/\-$get[0]\-/){
				#print "point $point\n";
					$point2=$point2."-".$get[1];
				}
				elsif($point=~/-$get[1]-/){
		                	$point2=$point2."-".$get[0];
		  		};
			$j++;
		};
	};
	
       	if($_=~/<POINTS>/){
       	
       	#print "point2 $point2\n";
       		if($point2 ne ''){
       		        $point2=~s/^-//;
       		        $point2=~s/-$//;
       			print DOC "> <POINTS>\n";
       			print DOC "$point2\n";
       			print DOC "\n";
       		};
       		$ecrit=0;
       		$pasvu=0;
       	}	
 	elsif($_=~/^\$\$\$\$/){
 		if($pasvu){
 		 	if($point2 ne ''){
 		 	        $point2=~s/^-//;
 		 	        $point2=~s/-$//;
       				print DOC "> <POINTS>\n";
       				print DOC "$point2\n";
       				print DOC "\n";
       			};
 		};
 		$nb++;
 		$point="";
 		$point2="";
 		$pasvu=1;
 		$comp=0;
 		$i=1;
 		$j=0; 	
 		print DOC $_;	
       	}
       	elsif($ecrit){
		print DOC $_;
	};
	$ecrit=1 if($get[0] eq '');
};
close(IN);
close(DOC);

#rename "tmp","$f";

print "$nb molecules dans $file\n";	

