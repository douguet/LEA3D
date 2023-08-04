#!/usr/bin/perl

	my($f)=@ARGV;
        if($f eq ''){
                print "usage: \n .pl <file operators.out>\n";
                exit;
        };

	$nbadd=0;
	$nbsup=0;
	$nbcross=0;
	$nbrep=0;
	$nbpermut=0;
	$i=0;
	open(MOL,"<$f");
        while(<MOL>){
		@get = split(' ',$_);
		$i++;
		$nbadd++ if($get[2] eq "add" && $get[3] > 0);
		$nbsup++ if($get[2] eq "suppress" && $get[3] > 0);
		$nbcross++ if($get[2] eq "crossover" && $get[3] > 0);
		$nbrep++ if($get[2] eq "replace" && $get[3] > 0);
		$nbpermut++ if($get[2] eq "permutation" && $get[3] > 0);
	};
	close(MOL);
	$pc=($nbcross/$i)*100;
	$pa=($nbadd/$i)*100;
	$ps=($nbsup/$i)*100;
	$pr=($nbrep/$i)*100;
	$pp=($nbpermut/$i)*100;
	print "crossover $pc add $pa supp $ps replace $pr permutation $pp\n";
