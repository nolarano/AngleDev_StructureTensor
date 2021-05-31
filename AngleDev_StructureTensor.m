%% Analyze vector fields find angle deviation and structure tensor, k
% Author: Sarah St. Pierre (last update 7//27/20)

% Download OrientationJ plug-in for ImageJ
% Run vector field for fiber images (make sure deg is chosen, not rad)
% Save results as .csv
% Use circle tool (hold for circle) and measure centroid location, radii
% length
% Put these measurements below last row in .csv file and save in .xlsx format

% Run section 1 for all 5 files; double-check Excel columns match variables
% Run section 2 and save meanrings and SEMrings
% Run section 4 (don't close histogram windows!)
    % >> ring# (1 through n) check bounds of bins (run in command window)
    % if the bounds are not [-96 90] you need to use alpharad# with
    % appropriate bounds for k for that ring#
    % if the bounds are not already in the list how to generate a new list:
        % for some bounds [#1 #2], take (|#1| + #2)/20 = ans1
        % ans1 / 2 = ans2
        % run (#1-ans2):ans1:(#2-ans2) to get new vector with alphadeg#
        % convert to alpharad#
% Run section 5 and check to make sure each k# has the correct alpharad#
%    corresponding to the bounds checked above; save k
        

%% Section 1
file1 = xlsread('tubulin5.xlsx'); % file name
dx = file1(1:12610,4); % dx from OrientationJ vector field results
dy = file1(1:12610,5); % dy "
coords = file1(1:12610,1:2); % x, y from vector field
xnot = file1(end,4)*72; % pattern centroid x, y location double-check scale is consistent between vector field results and circle measurement
ynot = file1(end,5)*72; % the 72 is scale correction for inch to pixel
ctrdist = sqrt((coords(:,1)-xnot).^2 + (coords(:,2)-ynot).^2);
cutoff = file1(end,7)*72/2; % boundary of pattern radius, check scale
outbound = find(ctrdist > cutoff);
ctrdist(outbound) = [];
dx(outbound) = [];
dy(outbound) = [];
coords(outbound,:) = [];
normdist = ctrdist/max(ctrdist);

% Angle Deviation
dotprod = (coords(:,1)-xnot).*dx + (coords(:,2)-ynot).*dy;
magvec = sqrt((coords(:,1)-xnot).^2 + (coords(:,2)-ynot).^2).*sqrt(dx.^2 + dy.^2);
angle = acos(dotprod./magvec);

% sort into rings
angdegrees = angle*180/pi;
histangle = angdegrees;
for i = 1:length(angdegrees)
    if angdegrees(i) > 90
        angdegrees(i) = 180 - angdegrees(i);
        histangle(i) = -angdegrees(i);
    end
end

k = 0:1/33:1; % 0:0.1:1 for 10 rings
ring_5 = cell(1,33); %change name for each file
hist_5 = cell(1,33); %change name
for i = 1:(length(k)-1)
    a = find(normdist < k(i+1) & normdist >= k(i));
    ring_5{i} = angdegrees(a); %change name
    hist_5{i} = histangle(a); % change name
end

k(1) = [];

%% Section 2: Compile all data together
alldata = cell(1,33);
allhist = cell(1,33);
for i = 1:33
    alldata{i} = [ring_1{i}; ring_2{i}; ring_3{i}; ring_4{i}; ring_5{i}];
    allhist{i} = [hist_1{i}; hist_2{i}; hist_3{i}; hist_4{i}; hist_5{i}];
end
% mean and SEM
meanrings = zeros(1,length(k));
SEMrings = zeros(1,length(k));
for j = 1:33
    meanrings(j) = mean(alldata{j}); % mean of fiber angle deviation - save this
    SEMrings(j) = std(alldata{j})/sqrt(length(alldata{j})); % Standard Error of Mean - save this
end
%% Section 3: Plot mean of all 5 images (Optional)
% I recommend using the saved meanrings and SEMrings results and plotting in
% Excel not here
errorbar(k,meanrings,SEMrings)
ylim([0 90]);
ylabel('Angle Deviation (degrees)')
xlabel('Normalized Distance from Center')
title('Control Circle Actin Angle Deviation')

%% Section 4: Histogram
nbins = 20;
figure(1)
ring1 = histogram(allhist{1},nbins);
A1 = sum(ring1.Values);
%ylim([0 0.1])
figure(2)
ring2 = histogram(allhist{2},nbins);
A2 = sum(ring2.Values);
%ylim([0 0.1])
figure(3)
ring3 = histogram(allhist{3},nbins);
A3 = sum(ring3.Values);
%ylim([0 0.1])
figure(4)
ring4 = histogram(allhist{4},nbins);
A4 = sum(ring4.Values);
%ylim([0 0.1])
figure(5)
ring5 = histogram(allhist{5},nbins);
A5 = sum(ring5.Values);
%ylim([0 0.1])
figure(6)
ring6 = histogram(allhist{6},nbins);
A6 = sum(ring6.Values);
%ylim([0 0.1])
figure(7)
ring7 = histogram(allhist{7},nbins);
A7 = sum(ring7.Values);
%ylim([0 0.1])
figure(8)
ring8 = histogram(allhist{8},nbins);
A8 = sum(ring8.Values);
%ylim([0 0.1])
figure(9)
ring9 = histogram(allhist{9},nbins);
A9 = sum(ring9.Values);
%ylim([0 0.1])
figure(10)
ring10 = histogram(allhist{10},nbins);
A10 = sum(ring10.Values);
%ylim([0 0.1])

figure(11)
ring11 = histogram(allhist{11},nbins);
A11 = sum(ring11.Values);
%ylim([0 0.1])
figure(12)
ring12 = histogram(allhist{12},nbins);
A12 = sum(ring12.Values);
%ylim([0 0.1])
figure(13)
ring13 = histogram(allhist{13},nbins);
A13 = sum(ring13.Values);
%ylim([0 0.1])
figure(14)
ring14 = histogram(allhist{14},nbins);
A14 = sum(ring14.Values);
%ylim([0 0.1])
figure(15)
ring15 = histogram(allhist{15},nbins);
A15 = sum(ring15.Values);
%ylim([0 0.1])
figure(16)
ring16 = histogram(allhist{16},nbins);
A16 = sum(ring16.Values);
%ylim([0 0.1])
figure(17)
ring17 = histogram(allhist{17},nbins);
A17 = sum(ring17.Values);
%ylim([0 0.1])
figure(18)
ring18 = histogram(allhist{18},nbins);
A18 = sum(ring18.Values);
%ylim([0 0.1])
figure(19)
ring19 = histogram(allhist{19},nbins);
A19 = sum(ring19.Values);
%ylim([0 0.1])
figure(20)
ring20 = histogram(allhist{20},nbins);
A20 = sum(ring20.Values);
%ylim([0 0.1])

figure(21)
ring21 = histogram(allhist{21},nbins);
A21 = sum(ring21.Values);
%ylim([0 0.1])
figure(22)
ring22 = histogram(allhist{22},nbins);
A22 = sum(ring22.Values);
%ylim([0 0.1])
figure(23)
ring23 = histogram(allhist{23},nbins);
A23 = sum(ring23.Values);
%ylim([0 0.1])
figure(24)
ring24 = histogram(allhist{24},nbins);
A24 = sum(ring24.Values);
%ylim([0 0.1])
figure(25)
ring25 = histogram(allhist{25},nbins);
A25 = sum(ring25.Values);
%ylim([0 0.1])
figure(26)
ring26 = histogram(allhist{26},nbins);
A26 = sum(ring26.Values);
%ylim([0 0.1])
figure(27)
ring27 = histogram(allhist{27},nbins);
A27 = sum(ring27.Values);
%ylim([0 0.1])
figure(28)
ring28 = histogram(allhist{28},nbins);
A28 = sum(ring28.Values);
%ylim([0 0.1])
figure(29)
ring29 = histogram(allhist{29},nbins);
A29 = sum(ring29.Values);
%ylim([0 0.1])
figure(30)
ring30 = histogram(allhist{30},nbins);
A30 = sum(ring30.Values);
%ylim([0 0.1])

figure(31)
ring31 = histogram(allhist{31},nbins);
A31 = sum(ring31.Values);
%ylim([0 0.1])
figure(32)
ring32 = histogram(allhist{32},nbins);
A32 = sum(ring32.Values);
% %ylim([0 0.1])
figure(33)
ring33 = histogram(allhist{33},nbins);
A33 = sum(ring33.Values);

%% Section 5: Finding K - Structure Tensor
alphadeg1 = [-87.0500000000000,-79.1500000000000,-71.2500000000000,-63.3500000000000,-55.4500000000000,-47.5500000000000,-39.6500000000000,-31.7500000000000,-23.8500000000000,-15.9500000000000,-8.05000000000000,-0.149999999999999,7.75000000000000,15.6500000000000,23.5500000000000,31.4500000000000,39.3500000000000,47.2500000000000,55.1500000000000,63.0500000000000];
alpharad1 = alphadeg1*pi/180;
alphadeg2 = [-93.7250000000000,-84.3855263157895,-75.0460526315790,-65.7065789473684,-56.3671052631579,-47.0276315789474,-37.6881578947368,-28.3486842105263,-19.0092105263158,-9.66973684210525,-0.330263157894748,9.00921052631578,18.3486842105263,27.6881578947368,37.0276315789474,46.3671052631579,55.7065789473684,65.0460526315790,74.3855263157895,83.7250000000000];
alpharad2 = alphadeg2*pi/180;
% [-96 84] - refers to below #
alphadeg3 = [-91.5000000000000,-82.5000000000000,-73.5000000000000,-64.5000000000000,-55.5000000000000,-46.5000000000000,-37.5000000000000,-28.5000000000000,-19.5000000000000,-10.5000000000000,-1.50000000000000,7.50000000000000,16.5000000000000,25.5000000000000,34.5000000000000,43.5000000000000,52.5000000000000,61.5000000000000,70.5000000000000,79.5000000000000];
alpharad3 = alphadeg3*pi/180;
% [-96 88]
alphadeg4 = [-91.4000000000000,-82.2000000000000,-73,-63.8000000000000,-54.6000000000000,-45.4000000000000,-36.2000000000000,-27,-17.8000000000000,-8.59999999999999,0.599999999999994,9.80000000000000,19,28.2000000000000,37.4000000000000,46.6000000000000,55.8000000000000,65,74.2000000000000,83.4000000000000];
alpharad4 = alphadeg4*pi/180;
alphadeg5 = [-83.8500000000000,-75.5500000000000,-67.2500000000000,-58.9500000000000,-50.6500000000000,-42.3500000000000,-34.0500000000000,-25.7500000000000,-17.4500000000000,-9.14999999999999,-0.850000000000009,7.44999999999999,15.7500000000000,24.0500000000000,32.3500000000000,40.6500000000000,48.9500000000000,57.2500000000000,65.5500000000000,73.8500000000000];
alpharad5 = alphadeg5*pi/180;
% [-84 62]
alphadeg6 = [-80.3500000000000,-73.0500000000000,-65.7500000000000,-58.4500000000000,-51.1500000000000,-43.8500000000000,-36.5500000000000,-29.2500000000000,-21.9500000000000,-14.6500000000000,-7.35000000000000,-0.0499999999999972,7.25000000000000,14.5500000000000,21.8500000000000,29.1500000000000,36.4500000000000,43.7500000000000,51.0500000000000,58.3500000000000];
alpharad6 = alphadeg6*pi/180;
% [-88 90]
alphadeg8 = [-83.5500000000000,-74.6500000000000,-65.7500000000000,-56.8500000000000,-47.9500000000000,-39.0500000000000,-30.1500000000000,-21.2500000000000,-12.3500000000000,-3.44999999999999,5.44999999999999,14.3500000000000,23.2500000000000,32.1500000000000,41.0500000000000,49.9500000000000,58.8500000000000,67.7500000000000,76.6500000000000,85.5500000000000];
alpharad8 = alphadeg8*pi/180;
% [-84 74]
alphadeg9 = [-80.0500000000000,-72.1500000000000,-64.2500000000000,-56.3500000000000,-48.4500000000000,-40.5500000000000,-32.6500000000000,-24.7500000000000,-16.8500000000000,-8.94999999999999,-1.05000000000001,6.84999999999999,14.7500000000000,22.6500000000000,30.5500000000000,38.4500000000000,46.3500000000000,54.2500000000000,62.1500000000000,70.0500000000000];
alpharad9 = alphadeg9*pi/180;
% [-88 82]
alphadeg10 = [-83.7500000000000,-75.2500000000000,-66.7500000000000,-58.2500000000000,-49.7500000000000,-41.2500000000000,-32.7500000000000,-24.2500000000000,-15.7500000000000,-7.25000000000000,1.25000000000000,9.75000000000000,18.2500000000000,26.7500000000000,35.2500000000000,43.7500000000000,52.2500000000000,60.7500000000000,69.2500000000000,77.7500000000000];
alpharad10 = alphadeg10*pi/180;
% [-88 84]
alphadeg11 = [-83.7000000000000,-75.1000000000000,-66.5000000000000,-57.9000000000000,-49.3000000000000,-40.7000000000000,-32.1000000000000,-23.5000000000000,-14.9000000000000,-6.30000000000001,2.30000000000001,10.9000000000000,19.5000000000000,28.1000000000000,36.7000000000000,45.3000000000000,53.9000000000000,62.5000000000000,71.1000000000000,79.7000000000000];
alpharad11 = alphadeg11*pi/180;
% [-66 62]
alphadeg12 = [-62.8000000000000,-56.4000000000000,-50,-43.6000000000000,-37.2000000000000,-30.8000000000000,-24.4000000000000,-18.0000000000000,-11.6000000000000,-5.20000000000000,1.20000000000000,7.59999999999999,14.0000000000000,20.4000000000000,26.8000000000000,33.2000000000000,39.6000000000000,46,52.4000000000000,58.8000000000000];
alpharad12 = alphadeg12*pi/180;
% [-88 86]
alphadeg13 = [-83.6500000000000,-74.9500000000000,-66.2500000000000,-57.5500000000000,-48.8500000000000,-40.1500000000000,-31.4500000000000,-22.7500000000000,-14.0500000000000,-5.35000000000001,3.35000000000001,12.0500000000000,20.7500000000000,29.4500000000000,38.1500000000000,46.8500000000000,55.5500000000000,64.2500000000000,72.9500000000000,81.6500000000000];
alpharad13 = alphadeg13*pi/180;
% [-96 88]
alphadeg14 = [-91.4500000000000,-82.3500000000000,-73.2500000000000,-64.1500000000000,-55.0500000000000,-45.9500000000000,-36.8500000000000,-27.7500000000000,-18.6500000000000,-9.55000000000001,-0.449999999999989,8.65000000000001,17.7500000000000,26.8500000000000,35.9500000000000,45.0500000000000,54.1500000000000,63.2500000000000,72.3500000000000,81.4500000000000];
alpharad14 = alphadeg14*pi/180;
% [-84 76]
alphadeg15 = [-80,-72,-64,-56,-48,-40,-32,-24,-16,-8,0,8,16,24,32,40,48,56,64,72];
alpharad15 = alphadeg15*pi/180;
% [-84 82]
alphadeg16 = [-79.8500000000000,-71.5500000000000,-63.2500000000000,-54.9500000000000,-46.6500000000000,-38.3500000000000,-30.0500000000000,-21.7500000000000,-13.4500000000000,-5.14999999999999,3.14999999999999,11.4500000000000,19.7500000000000,28.0500000000000,36.3500000000000,44.6500000000000,52.9500000000000,61.2500000000000,69.5500000000000,77.8500000000000];
alpharad16 = alphadeg16*pi/180;
% [-88 88]
alphadeg17 = [-83.6000000000000,-74.8000000000000,-66,-57.2000000000000,-48.4000000000000,-39.6000000000000,-30.8000000000000,-22.0000000000000,-13.2000000000000,-4.39999999999999,4.39999999999999,13.2000000000000,22.0000000000000,30.8000000000000,39.6000000000000,48.4000000000000,57.2000000000000,66,74.8000000000000,83.6000000000000];
alpharad17 = alphadeg17*pi/180;
% [-77 79]
alphadeg18 = [-73.1000000000000,-65.3000000000000,-57.5000000000000,-49.7000000000000,-41.9000000000000,-34.1000000000000,-26.3000000000000,-18.5000000000000,-10.7000000000000,-2.89999999999999,4.89999999999999,12.7000000000000,20.5000000000000,28.3000000000000,36.1000000000000,43.9000000000000,51.7000000000000,59.5000000000000,67.3000000000000,75.1000000000000];
alpharad18 = alphadeg18*pi/180;
% [-88 80]
alphadeg19 = [-83.8000000000000,-75.4000000000000,-67,-58.6000000000000,-50.2000000000000,-41.8000000000000,-33.4000000000000,-25.0000000000000,-16.6000000000000,-8.19999999999999,0.199999999999989,8.59999999999999,17.0000000000000,25.4000000000000,33.8000000000000,42.2000000000000,50.6000000000000,59,67.4000000000000,75.8000000000000];
alpharad19 = alphadeg19*pi/180;
% [-91 75]
alphadeg20 = [-86.8500000000000,-78.5500000000000,-70.2500000000000,-61.9500000000000,-53.6500000000000,-45.3500000000000,-37.0500000000000,-28.7500000000000,-20.4500000000000,-12.1500000000000,-3.85000000000001,4.44999999999999,12.7500000000000,21.0500000000000,29.3500000000000,37.6500000000000,45.9500000000000,54.2500000000000,62.5500000000000,70.8500000000000];
alpharad20 = alphadeg20*pi/180;
alphadeg21 = [-91.6000000000000,-82.8000000000000,-74,-65.2000000000000,-56.4000000000000,-47.6000000000000,-38.8000000000000,-30.0000000000000,-21.2000000000000,-12.4000000000000,-3.60000000000001,5.19999999999999,14.0000000000000,22.8000000000000,31.6000000000000,40.4000000000000,49.2000000000000,58.0000000000000,66.8000000000000,75.6000000000000];
alpharad21 = alphadeg21*pi/180;
% [-96 90]
alphadeg = [-91.35 -82.05 -72.75 -63.45 -54.15 -44.85 -35.55 -26.25 -16.95 -7.65 1.65 10.95 20.25 29.55 38.85 48.15 57.45 ...
    66.75 76.05 85.35];
alpharad = alphadeg*pi/180;
k = zeros(1,33);
for i = 1:20
    k(1) = k(1) + ring1(1).Values(i)*pi/A1*(sin(alpharad(i))).^2; % change alpharad# corresponding to ring# bounds
    k(2) = k(2) + ring2(1).Values(i)*pi/A2*(sin(alpharad(i))).^2;
    k(3) = k(3) + ring3(1).Values(i)*pi/A3*(sin(alpharad(i))).^2;
    k(4) = k(4) + ring4(1).Values(i)*pi/A4*(sin(alpharad(i))).^2;
    k(5) = k(5) + ring5(1).Values(i)*pi/A5*(sin(alpharad(i))).^2;
    k(6) = k(6) + ring6(1).Values(i)*pi/A6*(sin(alpharad(i))).^2;
    k(7) = k(7) + ring7(1).Values(i)*pi/A7*(sin(alpharad(i))).^2;
    k(8) = k(8) + ring8(1).Values(i)*pi/A8*(sin(alpharad(i))).^2;
    k(9) = k(9) + ring9(1).Values(i)*pi/A9*(sin(alpharad(i))).^2;
    k(10) = k(10) + ring10(1).Values(i)*pi/A10*(sin(alpharad(i))).^2;
    
    k(11) = k(11) + ring11(1).Values(i)*pi/A11*(sin(alpharad(i))).^2;
    k(12) = k(12) + ring12(1).Values(i)*pi/A12*(sin(alpharad(i))).^2;
    k(13) = k(13) + ring13(1).Values(i)*pi/A13*(sin(alpharad(i))).^2;
    k(14) = k(14) + ring14(1).Values(i)*pi/A14*(sin(alpharad(i))).^2;
    k(15) = k(15) + ring15(1).Values(i)*pi/A15*(sin(alpharad(i))).^2;
    k(16) = k(16) + ring16(1).Values(i)*pi/A16*(sin(alpharad(i))).^2;
    k(17) = k(17) + ring17(1).Values(i)*pi/A17*(sin(alpharad(i))).^2;
    k(18) = k(18) + ring18(1).Values(i)*pi/A18*(sin(alpharad(i))).^2;
    k(19) = k(19) + ring19(1).Values(i)*pi/A19*(sin(alpharad(i))).^2;
    k(20) = k(20) + ring20(1).Values(i)*pi/A20*(sin(alpharad(i))).^2;
    
    k(21) = k(21) + ring21(1).Values(i)*pi/A21*(sin(alpharad(i))).^2;
    k(22) = k(22) + ring22(1).Values(i)*pi/A22*(sin(alpharad(i))).^2;
    k(23) = k(23) + ring23(1).Values(i)*pi/A23*(sin(alpharad(i))).^2;
    k(24) = k(24) + ring24(1).Values(i)*pi/A24*(sin(alpharad(i))).^2;
    k(25) = k(25) + ring25(1).Values(i)*pi/A25*(sin(alpharad(i))).^2;
    k(26) = k(26) + ring26(1).Values(i)*pi/A26*(sin(alpharad(i))).^2;
    k(27) = k(27) + ring27(1).Values(i)*pi/A27*(sin(alpharad(i))).^2;
    k(28) = k(28) + ring28(1).Values(i)*pi/A28*(sin(alpharad(i))).^2;
    k(29) = k(29) + ring29(1).Values(i)*pi/A29*(sin(alpharad(i))).^2;
    k(30) = k(30) + ring30(1).Values(i)*pi/A30*(sin(alpharad(i))).^2;
    
    k(31) = k(31) + ring31(1).Values(i)*pi/A31*(sin(alpharad(i))).^2;
    k(32) = k(32) + ring32(1).Values(i)*pi/A32*(sin(alpharad(i))).^2;
    k(33) = k(33) + ring33(1).Values(i)*pi/A33*(sin(alpharad(i))).^2;
end

k_36hr_tubulin_33 = k*1/pi; % structure tensor - save this, make sure to rename
