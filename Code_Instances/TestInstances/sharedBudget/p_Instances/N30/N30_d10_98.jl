global arcs = [1 21; 2 30; 3 14; 3 24; 3 29; 4 6; 4 8; 5 3; 5 9; 5 26; 5 28; 6 23; 7 6; 7 20; 7 26; 7 28; 8 10; 8 12; 8 18; 8 25; 8 27; 9 3; 9 25; 10 8; 10 9; 10 13; 10 24; 11 8; 11 19; 12 24; 12 25; 12 26; 13 26; 14 13; 15 2; 15 4; 15 9; 15 17; 15 25; 15 27; 16 2; 16 25; 17 24; 17 29; 18 2; 18 3; 18 15; 18 22; 18 29; 19 20; 20 29; 21 7; 21 24; 21 29; 22 2; 22 5; 22 16; 23 7; 23 8; 23 11; 23 16; 23 21; 23 30; 24 21; 24 25; 24 28; 25 13; 25 14; 25 17; 26 25; 27 7; 27 14; 27 17; 27 25; 27 29; 28 5; 28 9; 28 23; 28 26; 29 9; 29 11]
global d_x = [6.0, 4.0, 5.0, 8.0, 1.0, 8.0, 8.0, 2.0, 9.0, 10.0, 6.0, 1.0, 4.0, 8.0, 9.0, 9.0, 3.0, 8.0, 6.0, 5.0, 10.0, 1.0, 6.0, 2.0, 7.0, 7.0, 9.0, 10.0, 4.0, 4.0, 7.0, 4.0, 7.0, 1.0, 6.0, 5.0, 10.0, 8.0, 6.0, 9.0, 5.0, 10.0, 6.0, 10.0, 1.0, 5.0, 9.0, 2.0, 5.0, 5.0, 10.0, 8.0, 8.0, 8.0, 8.0, 5.0, 10.0, 10.0, 3.0, 8.0, 6.0, 8.0, 9.0, 10.0, 8.0, 2.0, 10.0, 7.0, 9.0, 1.0, 7.0, 3.0, 8.0, 1.0, 6.0, 5.0, 8.0, 10.0, 5.0, 9.0, 6.0]
global b_x = 5
global d_y = [4.0, 3.0, 6.0, 1.0, 10.0, 7.0, 10.0, 9.0, 2.0, 7.0, 7.0, 2.0, 8.0, 2.0, 8.0, 9.0, 8.0, 10.0, 2.0, 10.0, 7.0, 7.0, 5.0, 1.0, 4.0, 4.0, 3.0, 8.0, 1.0, 2.0, 10.0, 9.0, 10.0, 7.0, 4.0, 8.0, 4.0, 8.0, 9.0, 3.0, 1.0, 1.0, 2.0, 1.0, 3.0, 3.0, 10.0, 6.0, 7.0, 8.0, 6.0, 9.0, 7.0, 3.0, 9.0, 1.0, 1.0, 10.0, 1.0, 9.0, 4.0, 2.0, 2.0, 7.0, 9.0, 10.0, 3.0, 3.0, 2.0, 9.0, 3.0, 1.0, 4.0, 1.0, 9.0, 5.0, 1.0, 10.0, 10.0, 2.0, 6.0]
global b_y = 10
global p = [0.562, 0.657, 0.18, 0.184, 0.909, 0.416, 0.027, 0.387, 0.091, 0.107, 0.392, 0.456, 0.414, 0.499, 0.289, 0.116, 0.963, 0.22, 0.261, 0.734, 0.673, 0.767, 0.69, 0.856, 0.18, 0.253, 0.694, 0.767, 0.345, 0.297, 0.737, 0.561, 0.643, 0.182, 0.383, 0.317, 0.64, 0.89, 0.16, 0.672, 0.514, 0.275, 0.199, 0.803, 0.576, 0.285, 0.222, 0.016, 0.2, 0.989, 0.682, 0.922, 0.942, 0.563, 0.614, 0.524, 0.478, 0.714, 0.332, 0.525, 0.144, 0.656, 0.05, 0.436, 0.081, 0.339, 0.436, 0.472, 0.085, 0.089, 0.703, 0.059, 0.02, 0.798, 0.397, 0.521, 0.143, 0.082, 0.799, 0.31, 0.901]
global q = [0.995, 0.896, 0.437, 0.48, 0.999, 0.712, 0.105, 0.608, 0.379, 0.923, 0.999, 0.689, 0.942, 0.889, 0.316, 0.66, 0.967, 0.288, 0.552, 0.926, 0.85, 0.933, 0.701, 0.974, 0.58, 0.541, 0.997, 0.829, 0.846, 0.302, 0.916, 0.844, 0.8, 0.204, 0.957, 0.98, 0.827, 0.917, 0.312, 0.925, 0.639, 0.367, 0.438, 0.959, 0.953, 0.822, 0.461, 0.189, 0.761, 0.995, 0.949, 0.95, 0.984, 0.693, 0.789, 0.789, 0.515, 0.961, 0.923, 0.724, 0.578, 0.826, 0.702, 0.653, 0.778, 0.68, 0.78, 0.571, 0.234, 0.735, 0.897, 0.336, 0.703, 0.987, 0.46, 0.741, 0.537, 0.305, 0.969, 0.566, 0.933]
global origin = 1
global destination = 30
