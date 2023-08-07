#!/usr/bin/perl


$f='';
my($f)=@ARGV;

if($f eq ''){
	print "usage: replacex <file.sdf>\n";
	exit;
};


$file="tmp.sdf";	
$blanc=" ";

open(DOC,">$file");
open(IN,"<$f");
while(<IN>){

        @get = split(' ',$_);

        if ($get[3] eq "X"){
        	$get[3]="H";
        	
        	$addh=sprintf"%10s%10s%10s%1s%1s%1s%3s%3s%3s\n",$get[0],$get[1],$get[2],$blanc,$get[3],$blanc,$get[4],$get[5],$get[6];
				 			
        	print DOC "$addh";
	}
	else{
		print DOC $_;
	};
	
};
close(IN);
close(DOC);

rename "$file","$f";






