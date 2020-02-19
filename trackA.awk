BEGIN{
    OFS="\t"
}
{   name = $1","$2
    split($3,a,"_");
    chr = a[1];
    pos = a[2]-1;
    str = a[3];
    type = a[4];
    color = $6>0 ? "255,140,0" : "139,0,139";
    label = "dPSI("name")="$6",q="$9;
    
    if(type=="D" && str=="+" || type=="A" && str=="-") {
	print chr, pos-1, pos+1, label, 1, str, pos, pos+1, color, 1, 2, 0
    }
    else {
	print chr, pos, pos+2, label, 1, str, pos, pos+1, color, 1, 2, 0
    }
}
