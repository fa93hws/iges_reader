% Auto Generated
function[objnodes, objlines, objsurfs, objregions, regioncolors,ishardPt] = spline_test_1_input()

objnodes = [
	-586.16847		,			 126.77609;		%	pt1
	-675.98915		,			 126.77609;		%	pt2
	-675.98915		,			 22.49699;		%	pt3
	-637.29958		,			 22.49699;		%	pt4
	-638.61852		,			 35.36356;		%	pt5
	-638.11214		,			 48.30480;		%	pt6
	-637.11546		,			 54.79928;		%	pt7
	-635.59164		,			 61.30743;		%	pt8
	-630.86818		,			 74.35815;		%	pt9
	-624.77962		,			 85.80653;		%	pt10
	-616.73336		,			 97.27265;		%	pt11
	-603.30798		,			 112.02653;		%	pt12
]; % 12pts in total

ishardPt = [1 1 1 1 0 0 0 0 0 0 0 0 ];
ishardPt = ishardPt==1;
objlines = cell(4,1) ;
objlines{1} 	= 	[1,2];
objlines{2} 	= 	[2,3];
objlines{3} 	= 	[3,4];
objlines{4} 	= 	[4,5,6,7,8,9,10,11,12,1];

objsurfs = cell(1,1) ;
objsurfs{1} 	 = 	[	1,2,3,4;
						1,1,1,1]';

objregions = cell(1,1) ;
objregions{1} 	 = 	[1];

regioncolors=[1];

%%			 regioncolors - layers
%color	layer_name
% 1		1H0
end
