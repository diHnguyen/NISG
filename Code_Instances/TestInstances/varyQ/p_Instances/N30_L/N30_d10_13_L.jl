global arcs = [1 12; 1 13; 1 22; 2 10; 2 14; 2 24; 2 26; 3 14; 3 16; 3 24; 3 26; 3 28; 4 28; 5 7; 5 9; 5 28; 6 9; 6 17; 6 24; 7 4; 7 12; 7 25; 8 11; 8 17; 9 15; 10 3; 11 9; 11 12; 11 13; 11 23; 11 30; 12 11; 12 20; 13 9; 13 10; 13 16; 13 23; 13 27; 13 30; 14 18; 15 8; 15 16; 15 28; 15 30; 16 7; 16 9; 16 12; 16 13; 16 22; 16 29; 17 2; 17 23; 18 5; 18 23; 19 8; 19 12; 19 14; 20 12; 21 6; 22 7; 22 15; 22 19; 23 21; 23 24; 23 25; 24 2; 24 7; 24 12; 25 10; 25 19; 26 15; 26 21; 26 22; 26 23; 27 5; 27 26; 28 8; 28 18; 28 27; 29 26; 29 27; 29 30]
global d_x = [10.0, 2.0, 10.0, 8.0, 4.0, 3.0, 1.0, 7.0, 10.0, 8.0, 1.0, 9.0, 5.0, 6.0, 3.0, 8.0, 7.0, 7.0, 6.0, 8.0, 9.0, 3.0, 7.0, 5.0, 4.0, 8.0, 4.0, 4.0, 1.0, 2.0, 4.0, 6.0, 2.0, 4.0, 9.0, 2.0, 8.0, 9.0, 2.0, 9.0, 3.0, 2.0, 4.0, 9.0, 10.0, 10.0, 9.0, 9.0, 7.0, 1.0, 8.0, 4.0, 3.0, 1.0, 10.0, 9.0, 4.0, 8.0, 7.0, 1.0, 4.0, 10.0, 1.0, 9.0, 5.0, 3.0, 8.0, 6.0, 9.0, 7.0, 10.0, 2.0, 8.0, 3.0, 1.0, 1.0, 5.0, 2.0, 9.0, 9.0, 3.0, 1.0]
global b_x = 5
global d_y = [5.0, 7.0, 5.0, 3.0, 3.0, 1.0, 3.0, 9.0, 8.0, 6.0, 3.0, 6.0, 5.0, 10.0, 6.0, 6.0, 6.0, 4.0, 8.0, 1.0, 3.0, 9.0, 7.0, 9.0, 2.0, 9.0, 9.0, 9.0, 6.0, 8.0, 4.0, 8.0, 5.0, 3.0, 3.0, 8.0, 10.0, 6.0, 6.0, 4.0, 1.0, 2.0, 4.0, 3.0, 7.0, 8.0, 7.0, 4.0, 6.0, 3.0, 6.0, 4.0, 1.0, 3.0, 8.0, 2.0, 5.0, 5.0, 7.0, 6.0, 4.0, 2.0, 2.0, 7.0, 7.0, 8.0, 6.0, 7.0, 6.0, 8.0, 4.0, 9.0, 7.0, 9.0, 6.0, 3.0, 6.0, 10.0, 3.0, 8.0, 10.0, 7.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.496, 0.499, 0.51, 0.509, 0.426, 0.431, 0.537, 0.566, 0.43, 0.513, 0.487, 0.58, 0.599, 0.44, 0.494, 0.495, 0.466, 0.461, 0.432, 0.556, 0.424, 0.416, 0.401, 0.402, 0.506, 0.57, 0.498, 0.584, 0.411, 0.411, 0.525, 0.517, 0.51, 0.491, 0.486, 0.577, 0.489, 0.423, 0.485, 0.596, 0.439, 0.438, 0.59, 0.557, 0.498, 0.517, 0.553, 0.485, 0.495, 0.416, 0.421, 0.535, 0.417, 0.486, 0.532, 0.562, 0.404, 0.521, 0.4, 0.413, 0.534, 0.429, 0.521, 0.418, 0.43, 0.408, 0.474, 0.565, 0.432, 0.542, 0.414, 0.424, 0.585, 0.425, 0.412, 0.561, 0.473, 0.575, 0.455, 0.558, 0.514, 0.406]
global origin = 1
global destination = 30