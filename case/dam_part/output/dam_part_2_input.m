% Auto Generated
function[objnodes, objlines, objsurfs, objregions, regioncolors,ishardPt] = dam_part_2_input()

objnodes = [
	-283.41310		,			-255.72487;		%	pt1
	-284.17430		,			-251.89804;		%	pt2
	-286.34203		,			-248.65380;		%	pt3
	-289.58626		,			-246.48607;		%	pt4
	-293.41310		,			-245.72487;		%	pt5
	-297.23993		,			-246.48607;		%	pt6
	-300.48417		,			-248.65380;		%	pt7
	-302.65189		,			-251.89804;		%	pt8
	-303.41310		,			-255.72487;		%	pt9
	-303.41310		,			-278.72487;		%	pt10
	-300.41310		,			-278.72487;		%	pt11
	-300.41310		,			-275.72487;		%	pt12
	-283.41310		,			-275.72487;		%	pt13
	-528.41310		,			-245.72487;		%	pt14
	-571.14690		,			-245.72487;		%	pt15
	-571.14690		,			-290.72487;		%	pt16
	-558.41310		,			-290.72487;		%	pt17
	-558.41310		,			-267.72487;		%	pt18
	-557.65189		,			-263.89804;		%	pt19
	-555.48417		,			-260.65380;		%	pt20
	-552.23993		,			-258.48607;		%	pt21
	-548.41310		,			-257.72487;		%	pt22
	-544.58626		,			-258.48607;		%	pt23
	-541.34203		,			-260.65380;		%	pt24
	-539.17430		,			-263.89804;		%	pt25
	-538.41310		,			-267.72487;		%	pt26
	-538.41310		,			-287.72487;		%	pt27
	-555.41310		,			-287.72487;		%	pt28
	-555.41310		,			-290.72487;		%	pt29
	-528.41310		,			-290.72487;		%	pt30
	 2025.34015		,			-111.72487;		%	pt31
	 6.31639		,			-111.72487;		%	pt32
	-83.18361		,			-290.72487;		%	pt33
	-308.41310		,			-290.72487;		%	pt34
	-961.41310		,			-290.72487;		%	pt35
	-953.91310		,			-305.72487;		%	pt36
	-73.91310		,			-305.72487;		%	pt37
	 15.58690		,			-126.72487;		%	pt38
	 2015.89990		,			-126.72487;		%	pt39
	 2066.70379		,			-66.72487;		%	pt40
	-21.49514		,			-66.72487;		%	pt41
	-110.99514		,			-245.72487;		%	pt42
	-273.41310		,			-245.72487;		%	pt43
	-273.41310		,			-283.72487;		%	pt44
	-313.41310		,			-283.72487;		%	pt45
	-313.41310		,			-245.72487;		%	pt46
	 2041.70323		,			-85.72487;		%	pt47
	 2066.70379		,			-85.72487;		%	pt48
	 2066.70379		,			 94.50785;		%	pt49
	 1809.58875		,			 11.19512;		%	pt50
	 593.58690		,			 17.27513;		%	pt51
	 593.58690		,			 23.27513;		%	pt52
	-816.89757		,			 23.27513;		%	pt53
	-816.89757		,			-210.72487;		%	pt54
	-771.41310		,			-210.72487;		%	pt55
	-771.41310		,			-245.72487;		%	pt56
	-313.41310		,			-240.72487;		%	pt57
	-273.41310		,			-240.72487;		%	pt58
	 2066.70379		,			 134.27513;		%	pt59
	 593.58690		,			 134.27513;		%	pt60
	 593.58690		,			 65.27513;		%	pt61
	-771.41310		,			-268.72487;		%	pt62
	-816.41310		,			-268.72487;		%	pt63
	-816.41310		,			-245.72487;		%	pt64
	-938.91310		,			-245.72487;		%	pt65
	-943.91310		,			-265.72487;		%	pt66
	-973.91310		,			-265.72487;		%	pt67
	-781.41310		,			-258.72487;		%	pt68
	-781.41310		,			-233.22487;		%	pt69
	-782.36460		,			-228.44133;		%	pt70
	-785.07426		,			-224.38603;		%	pt71
	-789.12955		,			-221.67638;		%	pt72
	-793.91310		,			-220.72487;		%	pt73
	-798.69664		,			-221.67638;		%	pt74
	-802.75193		,			-224.38603;		%	pt75
	-805.46159		,			-228.44133;		%	pt76
	-806.41310		,			-233.22487;		%	pt77
	-806.41310		,			-258.72487;		%	pt78
	-803.41310		,			-258.72487;		%	pt79
	-803.41310		,			-255.72487;		%	pt80
	-784.41310		,			-255.72487;		%	pt81
	-784.41310		,			-258.72487;		%	pt82
	-871.66310		,			 23.27513;		%	pt83
	-861.41310		,			 65.27513;		%	pt84
	-861.41310		,			 64.27513;		%	pt85
	 593.58690		,			 219.27513;		%	pt86
	 536.58690		,			 219.27513;		%	pt87
	 536.58690		,			 211.27513;		%	pt88
	 575.58690		,			 211.27513;		%	pt89
	 575.58690		,			 158.27513;		%	pt90
	 559.58690		,			 158.27513;		%	pt91
	 559.58690		,			 126.27513;		%	pt92
	 533.58690		,			 126.27513;		%	pt93
	 533.58690		,			 65.27513;		%	pt94
	 522.58690		,			 219.27513;		%	pt95
	 522.58690		,			 211.27513;		%	pt96
	 509.58690		,			 198.27513;		%	pt97
	 509.58690		,			 126.27513;		%	pt98
	 525.58690		,			 126.27513;		%	pt99
	 525.58690		,			 65.27513;		%	pt100
	-721.41310		,			 124.27513;		%	pt101
	-721.41310		,			 144.27513;		%	pt102
	-722.17430		,			 148.10196;		%	pt103
	-724.34203		,			 151.34620;		%	pt104
	-727.58626		,			 153.51393;		%	pt105
	-731.41310		,			 154.27513;		%	pt106
	-735.23993		,			 153.51393;		%	pt107
	-738.48417		,			 151.34620;		%	pt108
	-740.65189		,			 148.10196;		%	pt109
	-741.41310		,			 144.27513;		%	pt110
	-741.41310		,			 121.27513;		%	pt111
	-738.41310		,			 121.27513;		%	pt112
	-738.41310		,			 124.27513;		%	pt113
	-731.41310		,			 399.27513;		%	pt114
	-731.41310		,			 419.27513;		%	pt115
	-732.17430		,			 423.10196;		%	pt116
	-734.34203		,			 426.34620;		%	pt117
	-737.58626		,			 428.51393;		%	pt118
	-741.41310		,			 429.27513;		%	pt119
	-745.23993		,			 428.51393;		%	pt120
	-748.48417		,			 426.34620;		%	pt121
	-750.65189		,			 423.10196;		%	pt122
	-751.41310		,			 419.27513;		%	pt123
	-751.41310		,			 396.27513;		%	pt124
	-748.41310		,			 396.27513;		%	pt125
	-748.41310		,			 399.27513;		%	pt126
	-681.41310		,			 564.27513;		%	pt127
	-681.41310		,			 664.27513;		%	pt128
	-784.41310		,			 664.27513;		%	pt129
	-784.41310		,			 88.27513;		%	pt130
	-369.41310		,			 88.27513;		%	pt131
	-645.41153		,			 456.27304;		%	pt132
	-861.41310		,			 664.27513;		%	pt133
	 494.58690		,			 219.27513;		%	pt134
	 494.58690		,			 88.27513;		%	pt135
	-681.41310		,			 744.27513;		%	pt136
	-661.41310		,			 764.27513;		%	pt137
	-661.41310		,			 794.27513;		%	pt138
	-911.41310		,			 794.27513;		%	pt139
	-911.41310		,			 756.77513;		%	pt140
	-861.41310		,			 706.77513;		%	pt141
]; % 141pts in total

ishardPt = [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 1 1 1 1 1 1 1 1 1 ];
ishardPt = ishardPt==1;
objlines = cell(120,1) ;
objlines{1} 	= 	[1,2,3,4,5,6,7,8,9];
objlines{2} 	= 	[9,10];
objlines{3} 	= 	[10,11];
objlines{4} 	= 	[11,12];
objlines{5} 	= 	[12,13];
objlines{6} 	= 	[13,1];
objlines{7} 	= 	[15,16];
objlines{8} 	= 	[17,18];
objlines{9} 	= 	[18,19,20,21,22,23,24,25,26];
objlines{10} 	= 	[26,27];
objlines{11} 	= 	[27,28];
objlines{12} 	= 	[28,29];
objlines{13} 	= 	[30,14];
objlines{14} 	= 	[31,32];
objlines{15} 	= 	[32,33];
objlines{16} 	= 	[33,34];
objlines{17} 	= 	[34,30];
objlines{18} 	= 	[35,36];
objlines{19} 	= 	[36,37];
objlines{20} 	= 	[37,38];
objlines{21} 	= 	[38,39];
objlines{22} 	= 	[39,31];
objlines{23} 	= 	[40,41];
objlines{24} 	= 	[41,42];
objlines{25} 	= 	[42,43];
objlines{26} 	= 	[43,44];
objlines{27} 	= 	[44,45];
objlines{28} 	= 	[45,46];
objlines{29} 	= 	[31,47];
objlines{30} 	= 	[47,48];
objlines{31} 	= 	[48,40];
objlines{32} 	= 	[49,50];
objlines{33} 	= 	[50,51];
objlines{34} 	= 	[51,52];
objlines{35} 	= 	[53,54];
objlines{36} 	= 	[54,55];
objlines{37} 	= 	[56,15];
objlines{38} 	= 	[46,57];
objlines{39} 	= 	[57,58];
objlines{40} 	= 	[58,43];
objlines{41} 	= 	[40,49];
objlines{42} 	= 	[49,59];
objlines{43} 	= 	[59,60];
objlines{44} 	= 	[61,52];
objlines{45} 	= 	[62,63];
objlines{46} 	= 	[63,64];
objlines{47} 	= 	[64,65];
objlines{48} 	= 	[65,66];
objlines{49} 	= 	[66,67];
objlines{50} 	= 	[67,35];
objlines{51} 	= 	[68,69];
objlines{52} 	= 	[69,70,71,72,73,74,75,76,77];
objlines{53} 	= 	[77,78];
objlines{54} 	= 	[78,79];
objlines{55} 	= 	[79,80];
objlines{56} 	= 	[80,81];
objlines{57} 	= 	[81,82];
objlines{58} 	= 	[82,68];
objlines{59} 	= 	[83,65];
objlines{60} 	= 	[84,85];
objlines{61} 	= 	[85,83];
objlines{62} 	= 	[86,87];
objlines{63} 	= 	[87,88];
objlines{64} 	= 	[88,89];
objlines{65} 	= 	[89,90];
objlines{66} 	= 	[90,91];
objlines{67} 	= 	[91,92];
objlines{68} 	= 	[92,93];
objlines{69} 	= 	[93,94];
objlines{70} 	= 	[87,95];
objlines{71} 	= 	[95,96];
objlines{72} 	= 	[96,97];
objlines{73} 	= 	[97,98];
objlines{74} 	= 	[98,99];
objlines{75} 	= 	[99,100];
objlines{76} 	= 	[101,102];
objlines{77} 	= 	[102,103,104,105,106,107,108,109,110];
objlines{78} 	= 	[110,111];
objlines{79} 	= 	[111,112];
objlines{80} 	= 	[112,113];
objlines{81} 	= 	[113,101];
objlines{82} 	= 	[114,115];
objlines{83} 	= 	[115,116,117,118,119,120,121,122,123];
objlines{84} 	= 	[123,124];
objlines{85} 	= 	[124,125];
objlines{86} 	= 	[125,126];
objlines{87} 	= 	[126,114];
objlines{88} 	= 	[127,128];
objlines{89} 	= 	[129,130];
objlines{90} 	= 	[130,131];
objlines{91} 	= 	[131,132];
objlines{92} 	= 	[132,127];
objlines{93} 	= 	[133,84];
objlines{94} 	= 	[95,134];
objlines{95} 	= 	[134,135];
objlines{96} 	= 	[135,131];
objlines{97} 	= 	[136,137];
objlines{98} 	= 	[137,138];
objlines{99} 	= 	[138,139];
objlines{100} 	= 	[139,140];
objlines{101} 	= 	[140,141];
objlines{102} 	= 	[141,133];
objlines{103} 	= 	[128,136];
objlines{104} 	= 	[15,14];
objlines{105} 	= 	[14,46];
objlines{106} 	= 	[30,29];
objlines{107} 	= 	[29,17];
objlines{108} 	= 	[17,16];
objlines{109} 	= 	[16,35];
objlines{110} 	= 	[52,53];
objlines{111} 	= 	[53,83];
objlines{112} 	= 	[55,56];
objlines{113} 	= 	[56,62];
objlines{114} 	= 	[61,60];
objlines{115} 	= 	[60,86];
objlines{116} 	= 	[61,94];
objlines{117} 	= 	[94,100];
objlines{118} 	= 	[100,84];
objlines{119} 	= 	[128,129];
objlines{120} 	= 	[129,133];

objsurfs = cell(17,1) ;
objsurfs{1} 	 = 	[	1,2,3,4,5,6;
						1,1,1,1,1,1]';
objsurfs{2} 	 = 	[	104,7,108,8,9,10,11,12,106,13;
						-1,1,-1,1,1,1,1,1,-1,1]';
objsurfs{3} 	 = 	[	14,15,16,17,106,107,108,109,18,19,20,21,22;
						1,1,1,1,1,1,1,1,1,1,1,1,1]';
objsurfs{4} 	 = 	[	23,24,25,26,27,28,105,13,17,16,15,14,29,30,31;
						1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1,1,1]';
objsurfs{5} 	 = 	[	32,33,34,110,35,36,112,37,104,105,38,39,40,25,24,23,41;
						1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,1]';
objsurfs{6} 	 = 	[	42,43,114,44,34,33,32;
						1,1,-1,1,-1,-1,-1]';
objsurfs{7} 	 = 	[	37,113,45,46,47,48,49,50,109,7;
						-1,1,1,1,1,1,1,1,-1,-1]';
objsurfs{8} 	 = 	[	51,52,53,54,55,56,57,58;
						1,1,1,1,1,1,1,1]';
objsurfs{9} 	 = 	[	111,59,47,46,45,113,112,36,35;
						1,1,-1,-1,-1,-1,-1,-1,-1]';
objsurfs{10} 	 = 	[	116,117,118,60,61,111,110,44;
						1,1,1,1,1,-1,-1,-1]';
objsurfs{11} 	 = 	[	62,63,64,65,66,67,68,69,116,114,115;
						1,1,1,1,1,1,1,1,-1,1,1]';
objsurfs{12} 	 = 	[	64,63,70,71,72,73,74,75,117,69,68,67,66,65;
						-1,-1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1]';
objsurfs{13} 	 = 	[	76,77,78,79,80,81;
						1,1,1,1,1,1]';
objsurfs{14} 	 = 	[	82,83,84,85,86,87;
						1,1,1,1,1,1]';
objsurfs{15} 	 = 	[	88,119,89,90,91,92;
						1,1,1,1,1,1]';
objsurfs{16} 	 = 	[	120,93,118,75,74,73,72,71,94,95,96,90,89;
						1,1,-1,-1,-1,-1,-1,-1,1,1,1,-1,-1]';
objsurfs{17} 	 = 	[	97,98,99,100,101,102,120,119,103;
						1,1,1,1,1,1,-1,-1,1]';

objregions = cell(17,1) ;
objregions{1} 	 = 	[1];
objregions{2} 	 = 	[2];
objregions{3} 	 = 	[3];
objregions{4} 	 = 	[4];
objregions{5} 	 = 	[5];
objregions{6} 	 = 	[6];
objregions{7} 	 = 	[7];
objregions{8} 	 = 	[8];
objregions{9} 	 = 	[9];
objregions{10} 	 = 	[10];
objregions{11} 	 = 	[11];
objregions{12} 	 = 	[12];
objregions{13} 	 = 	[13];
objregions{14} 	 = 	[14];
objregions{15} 	 = 	[15];
objregions{16} 	 = 	[16];
objregions{17} 	 = 	[17];

regioncolors=[3,9,7,8,11,6,8,3,7,9,10,11,3,3,4,5,2];

%%			 regioncolors - layers
%color	layer_name
% 2		4Hmat0
% 3		7Hvoidmat
% 4		4Hmat1
% 5		4Hmat2
% 6		4Hmat8
% 7		4Hmat3
% 8		4Hmat4
% 9		4Hmat5
% 10		4Hmat7
% 11		4Hmat6
end