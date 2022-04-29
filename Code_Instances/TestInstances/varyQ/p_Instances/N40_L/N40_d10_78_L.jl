global arcs = [1 14; 1 19; 1 21; 1 22; 1 32; 2 4; 2 7; 2 12; 2 14; 2 38; 2 39; 3 2; 3 11; 3 23; 4 9; 4 13; 4 16; 4 36; 5 14; 5 23; 5 32; 5 33; 6 4; 6 18; 6 22; 6 27; 7 2; 7 18; 7 34; 7 38; 7 40; 8 5; 8 22; 8 35; 8 36; 9 19; 10 21; 11 2; 11 6; 11 13; 11 16; 11 21; 12 6; 12 10; 12 15; 12 17; 12 37; 12 40; 13 3; 13 8; 13 12; 13 29; 13 33; 14 8; 14 16; 14 17; 15 4; 15 16; 15 18; 15 26; 15 27; 15 32; 16 5; 16 10; 16 27; 16 35; 17 16; 17 33; 18 3; 18 20; 18 22; 18 24; 19 16; 20 6; 20 7; 20 14; 20 24; 20 36; 20 37; 21 5; 21 16; 21 20; 21 33; 21 39; 21 40; 22 8; 22 16; 22 25; 22 36; 23 8; 23 12; 23 36; 23 38; 24 30; 24 31; 24 40; 25 13; 25 24; 25 27; 26 4; 26 11; 26 16; 26 27; 26 36; 27 2; 27 13; 27 16; 27 37; 28 5; 28 12; 28 16; 29 4; 29 15; 29 28; 29 38; 30 6; 30 16; 30 27; 30 38; 31 22; 32 5; 32 13; 32 16; 32 17; 32 22; 33 7; 33 15; 33 39; 34 6; 34 23; 34 25; 35 7; 35 11; 35 26; 35 31; 35 32; 35 40; 36 33; 36 39; 37 4; 37 11; 37 39; 37 40; 38 7; 38 8; 38 35; 39 2; 39 7; 39 10; 39 19; 39 29]
global d_x = [8.0, 3.0, 10.0, 5.0, 3.0, 7.0, 8.0, 5.0, 10.0, 6.0, 5.0, 6.0, 6.0, 7.0, 6.0, 8.0, 4.0, 8.0, 6.0, 10.0, 10.0, 6.0, 6.0, 10.0, 9.0, 7.0, 7.0, 5.0, 4.0, 8.0, 2.0, 2.0, 10.0, 2.0, 2.0, 4.0, 5.0, 8.0, 5.0, 3.0, 10.0, 9.0, 8.0, 4.0, 5.0, 2.0, 7.0, 2.0, 8.0, 1.0, 1.0, 1.0, 5.0, 10.0, 8.0, 6.0, 5.0, 3.0, 6.0, 3.0, 5.0, 4.0, 2.0, 3.0, 9.0, 1.0, 4.0, 1.0, 4.0, 1.0, 5.0, 8.0, 10.0, 6.0, 9.0, 7.0, 7.0, 7.0, 4.0, 10.0, 1.0, 10.0, 6.0, 8.0, 8.0, 6.0, 10.0, 10.0, 1.0, 8.0, 2.0, 9.0, 10.0, 10.0, 9.0, 7.0, 10.0, 10.0, 10.0, 10.0, 4.0, 8.0, 6.0, 1.0, 6.0, 4.0, 1.0, 8.0, 7.0, 10.0, 9.0, 4.0, 9.0, 4.0, 3.0, 4.0, 10.0, 3.0, 9.0, 5.0, 5.0, 2.0, 5.0, 3.0, 2.0, 10.0, 7.0, 4.0, 5.0, 5.0, 4.0, 8.0, 3.0, 4.0, 1.0, 6.0, 3.0, 7.0, 1.0, 10.0, 10.0, 7.0, 5.0, 4.0, 4.0, 3.0, 6.0, 7.0, 7.0, 7.0, 6.0]
global b_x = 5
global d_y = [8.0, 3.0, 7.0, 6.0, 3.0, 9.0, 3.0, 8.0, 7.0, 3.0, 6.0, 1.0, 8.0, 10.0, 4.0, 9.0, 3.0, 3.0, 4.0, 1.0, 4.0, 9.0, 1.0, 8.0, 6.0, 8.0, 4.0, 5.0, 7.0, 4.0, 5.0, 10.0, 5.0, 2.0, 1.0, 4.0, 4.0, 7.0, 10.0, 3.0, 10.0, 3.0, 8.0, 8.0, 6.0, 5.0, 8.0, 9.0, 8.0, 9.0, 4.0, 2.0, 2.0, 1.0, 2.0, 3.0, 5.0, 2.0, 5.0, 9.0, 10.0, 6.0, 7.0, 5.0, 10.0, 4.0, 3.0, 4.0, 8.0, 4.0, 9.0, 2.0, 6.0, 3.0, 7.0, 1.0, 5.0, 4.0, 1.0, 4.0, 10.0, 1.0, 5.0, 7.0, 2.0, 8.0, 10.0, 8.0, 2.0, 3.0, 7.0, 6.0, 6.0, 7.0, 4.0, 8.0, 7.0, 7.0, 6.0, 10.0, 1.0, 7.0, 1.0, 9.0, 10.0, 7.0, 3.0, 2.0, 8.0, 8.0, 9.0, 2.0, 4.0, 10.0, 6.0, 4.0, 2.0, 9.0, 3.0, 4.0, 2.0, 4.0, 2.0, 7.0, 6.0, 6.0, 3.0, 6.0, 8.0, 1.0, 9.0, 6.0, 3.0, 2.0, 7.0, 8.0, 9.0, 10.0, 1.0, 5.0, 1.0, 2.0, 5.0, 9.0, 3.0, 9.0, 9.0, 8.0, 3.0, 6.0, 4.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.529, 0.565, 0.513, 0.55, 0.458, 0.569, 0.572, 0.449, 0.472, 0.566, 0.481, 0.563, 0.443, 0.574, 0.441, 0.485, 0.502, 0.574, 0.491, 0.477, 0.575, 0.577, 0.545, 0.411, 0.436, 0.434, 0.428, 0.586, 0.569, 0.521, 0.561, 0.542, 0.568, 0.523, 0.571, 0.426, 0.475, 0.515, 0.559, 0.546, 0.525, 0.478, 0.482, 0.574, 0.444, 0.569, 0.54, 0.429, 0.432, 0.533, 0.577, 0.597, 0.502, 0.4, 0.516, 0.525, 0.582, 0.468, 0.554, 0.414, 0.488, 0.492, 0.475, 0.512, 0.426, 0.575, 0.435, 0.556, 0.453, 0.529, 0.594, 0.502, 0.425, 0.543, 0.487, 0.504, 0.489, 0.504, 0.531, 0.426, 0.56, 0.552, 0.418, 0.51, 0.503, 0.516, 0.577, 0.511, 0.533, 0.481, 0.426, 0.431, 0.429, 0.543, 0.436, 0.525, 0.566, 0.492, 0.451, 0.552, 0.42, 0.592, 0.533, 0.513, 0.509, 0.411, 0.586, 0.49, 0.455, 0.45, 0.433, 0.529, 0.522, 0.454, 0.551, 0.408, 0.415, 0.519, 0.456, 0.525, 0.495, 0.415, 0.507, 0.511, 0.591, 0.531, 0.572, 0.456, 0.408, 0.568, 0.578, 0.446, 0.42, 0.429, 0.494, 0.596, 0.565, 0.428, 0.449, 0.55, 0.598, 0.515, 0.597, 0.57, 0.559, 0.491, 0.506, 0.472, 0.495, 0.597, 0.439]
global origin = 1
global destination = 40