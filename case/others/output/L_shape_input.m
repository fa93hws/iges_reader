% Auto Generated
function[objnodes, objlines, objsurfs, objregions, regioncolors,ishardPt] = L_shape_input()

objnodes = [
	 190.10535		,			 151.78734;		%	pt1
	-161.49669		,			 151.78734;		%	pt2
	-161.49669		,			-139.49616;		%	pt3
	-73.21566		,			-139.49616;		%	pt4
	-73.21566		,			 43.79189;		%	pt5
	 190.10535		,			 43.79189;		%	pt6
]; % 6pts in total

ishardPt = [1 1 1 1 1 1 ];
ishardPt = ishardPt==1;
objlines = cell(6,1) ;
objlines{1} 	= 	[1,2];
objlines{2} 	= 	[2,3];
objlines{3} 	= 	[3,4];
objlines{4} 	= 	[4,5];
objlines{5} 	= 	[5,6];
objlines{6} 	= 	[6,1];

objsurfs = cell(1,1) ;
objsurfs{1} 	 = 	[	1,2,3,4,5,6;
						1,1,1,1,1,1]';

objregions = cell(1,1) ;
objregions{1} 	 = 	[1];

regioncolors=[1];

%%			 regioncolors - layers
%color	layer_name
% 1		1H0
end
