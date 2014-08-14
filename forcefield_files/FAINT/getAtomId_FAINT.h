int getAtomId_FAINT(string, string);

int getAtomId_FAINT(string rnam, string anam){
if(rnam.compare("ASP")==0){
	 if(anam.compare("OD2")==0) return 31;
	 if(anam.compare("O")==0) return 27;
	 if(anam.compare("OD1")==0) return 30;
	 if(anam.compare("C")==0) return 26;
	 if(anam.compare("N")==0) return 24;
	 if(anam.compare("CG")==0) return 29;
	 if(anam.compare("CB")==0) return 28;
	 if(anam.compare("CA")==0) return 25;
}

if(rnam.compare("PRO")==0){
	 if(anam.compare("CD")==0) return 120;
	 if(anam.compare("O")==0) return 117;
	 if(anam.compare("C")==0) return 116;
	 if(anam.compare("N")==0) return 114;
	 if(anam.compare("CG")==0) return 119;
	 if(anam.compare("CB")==0) return 118;
	 if(anam.compare("CA")==0) return 115;
}

if(rnam.compare("LYS")==0){
	 if(anam.compare("O")==0) return 89;
	 if(anam.compare("N")==0) return 86;
	 if(anam.compare("CB")==0) return 90;
	 if(anam.compare("CE")==0) return 93;
	 if(anam.compare("CD")==0) return 92;
	 if(anam.compare("C")==0) return 88;
	 if(anam.compare("CG")==0) return 91;
	 if(anam.compare("CA")==0) return 87;
	 if(anam.compare("NZ")==0) return 94;
}

if(rnam.compare("ILE")==0){
	 if(anam.compare("O")==0) return 73;
	 if(anam.compare("N")==0) return 70;
	 if(anam.compare("CB")==0) return 74;
	 if(anam.compare("CG2")==0) return 76;
	 if(anam.compare("CD1")==0) return 77;
	 if(anam.compare("C")==0) return 72;
	 if(anam.compare("CG1")==0) return 75;
	 if(anam.compare("CA")==0) return 71;
}

if(rnam.compare("TRP")==0){
	 if(anam.compare("O")==0) return 137;
	 if(anam.compare("CZ2")==0) return 145;
	 if(anam.compare("N")==0) return 134;
	 if(anam.compare("CB")==0) return 138;
	 if(anam.compare("CE2")==0) return 143;
	 if(anam.compare("NE1")==0) return 142;
	 if(anam.compare("CD1")==0) return 140;
	 if(anam.compare("CD2")==0) return 141;
	 if(anam.compare("CH2")==0) return 147;
	 if(anam.compare("C")==0) return 136;
	 if(anam.compare("CE3")==0) return 144;
	 if(anam.compare("CG")==0) return 139;
	 if(anam.compare("CZ3")==0) return 146;
	 if(anam.compare("CA")==0) return 135;
}

if(rnam.compare("CYS")==0){
	 if(anam.compare("O")==0) return 35;
	 if(anam.compare("C")==0) return 34;
	 if(anam.compare("N")==0) return 32;
	 if(anam.compare("CB")==0) return 36;
	 if(anam.compare("CA")==0) return 33;
	 if(anam.compare("SG")==0) return 37;
}

if(rnam.compare("GLY")==0){
	 if(anam.compare("O")==0) return 59;
	 if(anam.compare("C")==0) return 58;
	 if(anam.compare("N")==0) return 56;
	 if(anam.compare("CA")==0) return 57;
}

if(rnam.compare("PHE")==0){
	 if(anam.compare("O")==0) return 106;
	 if(anam.compare("N")==0) return 103;
	 if(anam.compare("CB")==0) return 107;
	 if(anam.compare("CE2")==0) return 112;
	 if(anam.compare("CD1")==0) return 109;
	 if(anam.compare("CD2")==0) return 110;
	 if(anam.compare("CZ")==0) return 113;
	 if(anam.compare("CE1")==0) return 111;
	 if(anam.compare("C")==0) return 105;
	 if(anam.compare("CG")==0) return 108;
	 if(anam.compare("CA")==0) return 104;
}

if(rnam.compare("GLN")==0){
	 if(anam.compare("O")==0) return 41;
	 if(anam.compare("NE2")==0) return 46;
	 if(anam.compare("N")==0) return 38;
	 if(anam.compare("OE1")==0) return 45;
	 if(anam.compare("CB")==0) return 42;
	 if(anam.compare("CD")==0) return 44;
	 if(anam.compare("C")==0) return 40;
	 if(anam.compare("CG")==0) return 43;
	 if(anam.compare("CA")==0) return 39;
}

if(rnam.compare("SER")==0){
	 if(anam.compare("O")==0) return 124;
	 if(anam.compare("OG")==0) return 126;
	 if(anam.compare("C")==0) return 123;
	 if(anam.compare("N")==0) return 121;
	 if(anam.compare("CB")==0) return 125;
	 if(anam.compare("CA")==0) return 122;
}

if(rnam.compare("ASN")==0){
	 if(anam.compare("O")==0) return 19;
	 if(anam.compare("N")==0) return 16;
	 if(anam.compare("CB")==0) return 20;
	 if(anam.compare("ND2")==0) return 23;
	 if(anam.compare("C")==0) return 18;
	 if(anam.compare("OD1")==0) return 22;
	 if(anam.compare("CG")==0) return 21;
	 if(anam.compare("CA")==0) return 17;
}

if(rnam.compare("PRB")==0){
	 if(anam.compare("DON")==0) return 169;
	 if(anam.compare("ACC")==0) return 170;
	 if(anam.compare("HYD")==0) return 167;
	 if(anam.compare("NEG")==0) return 172;
	 if(anam.compare("ARM")==0) return 168;
	 if(anam.compare("POS")==0) return 171;
}

if(rnam.compare("VAL")==0){
	 if(anam.compare("CG2")==0) return 166;
	 if(anam.compare("O")==0) return 163;
	 if(anam.compare("C")==0) return 162;
	 if(anam.compare("N")==0) return 160;
	 if(anam.compare("CG1")==0) return 165;
	 if(anam.compare("CB")==0) return 164;
	 if(anam.compare("CA")==0) return 161;
}

if(rnam.compare("LEU")==0){
	 if(anam.compare("O")==0) return 81;
	 if(anam.compare("N")==0) return 78;
	 if(anam.compare("CB")==0) return 82;
	 if(anam.compare("CD1")==0) return 84;
	 if(anam.compare("CD2")==0) return 85;
	 if(anam.compare("C")==0) return 80;
	 if(anam.compare("CG")==0) return 83;
	 if(anam.compare("CA")==0) return 79;
}

if(rnam.compare("TYR")==0){
	 if(anam.compare("O")==0) return 151;
	 if(anam.compare("N")==0) return 148;
	 if(anam.compare("CB")==0) return 152;
	 if(anam.compare("CE2")==0) return 157;
	 if(anam.compare("CD1")==0) return 154;
	 if(anam.compare("CD2")==0) return 155;
	 if(anam.compare("CZ")==0) return 158;
	 if(anam.compare("CE1")==0) return 156;
	 if(anam.compare("C")==0) return 150;
	 if(anam.compare("CG")==0) return 153;
	 if(anam.compare("OH")==0) return 159;
	 if(anam.compare("CA")==0) return 149;
}

if(rnam.compare("GLU")==0){
	 if(anam.compare("O")==0) return 50;
	 if(anam.compare("N")==0) return 47;
	 if(anam.compare("OE1")==0) return 54;
	 if(anam.compare("CB")==0) return 51;
	 if(anam.compare("OE2")==0) return 55;
	 if(anam.compare("CD")==0) return 53;
	 if(anam.compare("C")==0) return 49;
	 if(anam.compare("CG")==0) return 52;
	 if(anam.compare("CA")==0) return 48;
}

if(rnam.compare("ARG")==0){
	 if(anam.compare("O")==0) return 8;
	 if(anam.compare("N")==0) return 5;
	 if(anam.compare("CB")==0) return 9;
	 if(anam.compare("CD")==0) return 11;
	 if(anam.compare("CZ")==0) return 13;
	 if(anam.compare("C")==0) return 7;
	 if(anam.compare("CG")==0) return 10;
	 if(anam.compare("NE")==0) return 12;
	 if(anam.compare("NH2")==0) return 15;
	 if(anam.compare("NH1")==0) return 14;
	 if(anam.compare("CA")==0) return 6;
}

if(rnam.compare("THR")==0){
	 if(anam.compare("CG2")==0) return 133;
	 if(anam.compare("OG1")==0) return 132;
	 if(anam.compare("O")==0) return 130;
	 if(anam.compare("C")==0) return 129;
	 if(anam.compare("N")==0) return 127;
	 if(anam.compare("CB")==0) return 131;
	 if(anam.compare("CA")==0) return 128;
}

if(rnam.compare("ALA")==0){
	 if(anam.compare("O")==0) return 3;
	 if(anam.compare("C")==0) return 2;
	 if(anam.compare("N")==0) return 0;
	 if(anam.compare("CB")==0) return 4;
	 if(anam.compare("CA")==0) return 1;
}

if(rnam.compare("MET")==0){
	 if(anam.compare("O")==0) return 98;
	 if(anam.compare("N")==0) return 95;
	 if(anam.compare("CB")==0) return 99;
	 if(anam.compare("CE")==0) return 102;
	 if(anam.compare("C")==0) return 97;
	 if(anam.compare("CG")==0) return 100;
	 if(anam.compare("SD")==0) return 101;
	 if(anam.compare("CA")==0) return 96;
}

if(rnam.compare("HIS")==0){
	 if(anam.compare("O")==0) return 63;
	 if(anam.compare("NE2")==0) return 69;
	 if(anam.compare("N")==0) return 60;
	 if(anam.compare("CB")==0) return 64;
	 if(anam.compare("ND1")==0) return 66;
	 if(anam.compare("CD2")==0) return 67;
	 if(anam.compare("C")==0) return 62;
	 if(anam.compare("CE1")==0) return 68;
	 if(anam.compare("CG")==0) return 65;
	 if(anam.compare("CA")==0) return 61;
}


return(-1);}
