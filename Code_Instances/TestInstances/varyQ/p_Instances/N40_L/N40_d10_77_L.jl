global arcs = [1 7; 1 9; 1 11; 1 36; 1 37; 1 39; 2 18; 2 21; 2 23; 2 40; 3 6; 3 20; 3 21; 3 26; 3 30; 4 11; 4 16; 4 26; 4 30; 4 39; 5 10; 5 24; 6 3; 6 7; 6 34; 7 12; 7 18; 7 25; 8 5; 8 6; 8 18; 8 20; 8 33; 8 35; 8 37; 9 28; 10 14; 10 21; 10 24; 10 27; 11 7; 11 27; 11 33; 11 36; 11 38; 11 39; 12 2; 12 29; 13 9; 13 17; 13 24; 13 28; 13 39; 14 4; 14 8; 14 19; 14 22; 14 30; 14 37; 15 23; 15 25; 15 29; 16 20; 16 21; 16 27; 16 31; 16 38; 17 4; 17 7; 17 9; 17 12; 17 14; 17 29; 17 33; 18 13; 18 36; 18 40; 19 16; 19 23; 19 34; 19 35; 20 13; 20 19; 20 34; 20 39; 21 6; 21 20; 21 22; 21 31; 21 33; 22 8; 22 11; 22 32; 23 12; 23 36; 23 40; 24 12; 24 23; 24 30; 25 29; 25 39; 26 9; 26 23; 27 4; 27 11; 27 16; 27 22; 28 29; 29 36; 29 38; 30 12; 30 33; 30 40; 31 4; 31 8; 31 23; 31 35; 31 39; 32 10; 32 16; 32 18; 32 33; 33 19; 33 37; 34 2; 34 3; 34 7; 34 12; 35 3; 35 15; 35 19; 36 2; 36 3; 36 10; 36 16; 36 21; 37 4; 37 26; 38 11; 39 11; 39 12; 39 26; 39 29; 39 33]
global d_x = [6.0, 8.0, 7.0, 10.0, 9.0, 4.0, 9.0, 6.0, 6.0, 10.0, 1.0, 1.0, 10.0, 7.0, 7.0, 10.0, 3.0, 8.0, 3.0, 7.0, 4.0, 1.0, 7.0, 10.0, 7.0, 10.0, 10.0, 3.0, 9.0, 5.0, 7.0, 4.0, 8.0, 4.0, 3.0, 2.0, 1.0, 8.0, 4.0, 2.0, 8.0, 1.0, 8.0, 2.0, 2.0, 5.0, 10.0, 7.0, 4.0, 1.0, 4.0, 10.0, 3.0, 5.0, 8.0, 2.0, 4.0, 4.0, 2.0, 10.0, 8.0, 6.0, 8.0, 7.0, 1.0, 8.0, 8.0, 4.0, 5.0, 5.0, 4.0, 2.0, 9.0, 10.0, 9.0, 9.0, 5.0, 4.0, 5.0, 8.0, 9.0, 4.0, 3.0, 5.0, 2.0, 7.0, 8.0, 2.0, 9.0, 3.0, 9.0, 4.0, 6.0, 7.0, 8.0, 5.0, 3.0, 1.0, 9.0, 2.0, 2.0, 9.0, 10.0, 1.0, 10.0, 7.0, 4.0, 2.0, 8.0, 1.0, 1.0, 6.0, 2.0, 1.0, 4.0, 10.0, 4.0, 8.0, 8.0, 10.0, 6.0, 4.0, 9.0, 3.0, 2.0, 8.0, 2.0, 9.0, 7.0, 9.0, 1.0, 2.0, 3.0, 6.0, 7.0, 8.0, 2.0, 3.0, 10.0, 8.0, 5.0, 7.0, 3.0, 1.0]
global b_x = 5
global d_y = [10.0, 3.0, 3.0, 1.0, 9.0, 7.0, 1.0, 8.0, 10.0, 10.0, 6.0, 7.0, 4.0, 10.0, 6.0, 1.0, 7.0, 4.0, 8.0, 2.0, 1.0, 9.0, 8.0, 7.0, 10.0, 3.0, 9.0, 1.0, 8.0, 1.0, 8.0, 4.0, 6.0, 10.0, 10.0, 10.0, 2.0, 7.0, 5.0, 7.0, 7.0, 8.0, 5.0, 3.0, 9.0, 10.0, 7.0, 2.0, 9.0, 2.0, 4.0, 4.0, 2.0, 3.0, 5.0, 10.0, 5.0, 10.0, 3.0, 2.0, 3.0, 1.0, 5.0, 5.0, 9.0, 5.0, 4.0, 2.0, 5.0, 1.0, 8.0, 4.0, 8.0, 4.0, 2.0, 3.0, 7.0, 5.0, 9.0, 8.0, 6.0, 4.0, 4.0, 3.0, 4.0, 7.0, 4.0, 5.0, 4.0, 7.0, 3.0, 5.0, 10.0, 5.0, 1.0, 3.0, 3.0, 6.0, 4.0, 2.0, 10.0, 10.0, 6.0, 9.0, 5.0, 10.0, 9.0, 2.0, 8.0, 9.0, 6.0, 9.0, 2.0, 2.0, 4.0, 3.0, 10.0, 7.0, 4.0, 2.0, 6.0, 7.0, 8.0, 6.0, 1.0, 2.0, 10.0, 2.0, 5.0, 5.0, 10.0, 1.0, 8.0, 1.0, 1.0, 8.0, 10.0, 7.0, 2.0, 6.0, 6.0, 5.0, 2.0, 2.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.539, 0.545, 0.425, 0.508, 0.472, 0.469, 0.555, 0.488, 0.56, 0.443, 0.57, 0.54, 0.418, 0.432, 0.524, 0.52, 0.526, 0.57, 0.433, 0.428, 0.547, 0.42, 0.571, 0.544, 0.53, 0.478, 0.592, 0.452, 0.554, 0.451, 0.547, 0.459, 0.572, 0.402, 0.517, 0.483, 0.484, 0.41, 0.557, 0.442, 0.513, 0.526, 0.499, 0.483, 0.569, 0.496, 0.518, 0.534, 0.441, 0.446, 0.499, 0.567, 0.454, 0.451, 0.433, 0.468, 0.59, 0.429, 0.465, 0.45, 0.499, 0.58, 0.467, 0.592, 0.484, 0.494, 0.433, 0.577, 0.458, 0.532, 0.551, 0.586, 0.472, 0.508, 0.427, 0.586, 0.587, 0.522, 0.585, 0.526, 0.425, 0.448, 0.477, 0.516, 0.402, 0.46, 0.524, 0.516, 0.564, 0.401, 0.594, 0.447, 0.462, 0.497, 0.48, 0.52, 0.569, 0.535, 0.504, 0.456, 0.499, 0.411, 0.481, 0.411, 0.504, 0.557, 0.44, 0.457, 0.599, 0.502, 0.521, 0.508, 0.428, 0.449, 0.523, 0.544, 0.424, 0.564, 0.557, 0.414, 0.46, 0.587, 0.56, 0.4, 0.561, 0.402, 0.454, 0.508, 0.417, 0.498, 0.439, 0.557, 0.471, 0.435, 0.584, 0.581, 0.597, 0.467, 0.571, 0.426, 0.435, 0.414, 0.594, 0.531]
global origin = 1
global destination = 40