% Auto Generated
function[objnodes, objlines, objsurfs, objregions, regioncolors,ishardPt] = dam_part_surface_1_input()

objnodes = [
	 496.48094		,			-128.99956;		%	pt1
	 401.48094		,			-128.99956;		%	pt2
	 401.48094		,			-110.55802;		%	pt3
	 371.07378		,			-110.55802;		%	pt4
	 353.73094		,			-89.99956;		%	pt5
	 308.48094		,			-29.66623;		%	pt6
	 308.48094		,			 0.33377;		%	pt7
	 185.98094		,			 0.33377;		%	pt8
	 173.48094		,			 17.00044;		%	pt9
	-171.51906		,			 17.00044;		%	pt10
	-171.51906		,			-35.99956;		%	pt11
	-211.51906		,			-35.99956;		%	pt12
	-211.51906		,			-357.99956;		%	pt13
	-171.51906		,			-357.99956;		%	pt14
	-171.51906		,			-392.99956;		%	pt15
	 31.48094		,			-392.99956;		%	pt16
	 71.48094		,			-392.99956;		%	pt17
	 286.48094		,			-392.99956;		%	pt18
	 326.48094		,			-392.99956;		%	pt19
	 536.48094		,			-392.99956;		%	pt20
	 536.48094		,			-387.99956;		%	pt21
	 576.48094		,			-387.99956;		%	pt22
	 576.48094		,			-392.99956;		%	pt23
	 658.48094		,			-392.99956;		%	pt24
	 658.48094		,			-178.99956;		%	pt25
	 516.48094		,			-178.99956;		%	pt26
	 496.48094		,			-158.99956;		%	pt27
]; % 27pts in total

ishardPt = [1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1 1 ];
ishardPt = ishardPt==1;
objlines = cell(27,1) ;
objlines{1} 	= 	[1,2];
objlines{2} 	= 	[2,3];
objlines{3} 	= 	[3,4];
objlines{4} 	= 	[4,5];
objlines{5} 	= 	[5,6];
objlines{6} 	= 	[6,7];
objlines{7} 	= 	[7,8];
objlines{8} 	= 	[8,9];
objlines{9} 	= 	[9,10];
objlines{10} 	= 	[10,11];
objlines{11} 	= 	[11,12];
objlines{12} 	= 	[12,13];
objlines{13} 	= 	[13,14];
objlines{14} 	= 	[14,15];
objlines{15} 	= 	[15,16];
objlines{16} 	= 	[16,17];
objlines{17} 	= 	[17,18];
objlines{18} 	= 	[18,19];
objlines{19} 	= 	[19,20];
objlines{20} 	= 	[20,21];
objlines{21} 	= 	[21,22];
objlines{22} 	= 	[22,23];
objlines{23} 	= 	[23,24];
objlines{24} 	= 	[24,25];
objlines{25} 	= 	[25,26];
objlines{26} 	= 	[26,27];
objlines{27} 	= 	[27,1];

objsurfs = cell(1,1) ;
objsurfs{1} 	 = 	[	1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27;
						1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]';

objregions = cell(1,1) ;
objregions{1} 	 = 	[1];

regioncolors=[2];

%%			 regioncolors - layers
%color	layer_name
% 2		4Hmat2
end
