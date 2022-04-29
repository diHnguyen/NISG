global arcs = [1 2; 1 3; 1 10; 1 35; 2 5; 2 6; 2 22; 2 28; 3 8; 3 10; 3 15; 3 23; 3 24; 3 25; 4 5; 4 14; 4 28; 5 14; 5 30; 6 11; 7 4; 7 8; 7 21; 7 28; 8 3; 8 5; 8 15; 8 25; 8 29; 8 35; 9 3; 9 8; 9 11; 9 25; 10 8; 10 25; 10 28; 10 31; 11 9; 11 12; 11 20; 12 11; 12 17; 12 23; 12 29; 13 3; 13 7; 13 10; 13 28; 14 11; 14 28; 15 5; 15 11; 15 13; 15 22; 15 27; 15 31; 15 33; 16 4; 16 13; 16 18; 17 14; 17 22; 17 26; 17 33; 18 2; 18 7; 18 19; 18 20; 18 25; 18 34; 19 11; 19 12; 19 14; 19 26; 20 10; 20 13; 20 22; 20 23; 21 6; 21 10; 21 17; 21 20; 21 33; 22 12; 22 34; 23 2; 23 12; 23 24; 23 27; 24 16; 24 17; 24 29; 24 30; 24 31; 25 6; 25 29; 25 35; 26 2; 26 22; 26 29; 26 32; 27 21; 27 23; 27 34; 27 35; 28 6; 28 17; 28 30; 29 9; 30 19; 31 3; 31 10; 31 21; 32 4; 32 34; 33 31; 34 4; 34 5; 34 15; 34 16; 34 27; 34 29; 34 31]
global d_x = [5.0, 5.0, 9.0, 1.0, 8.0, 6.0, 6.0, 5.0, 4.0, 6.0, 6.0, 7.0, 7.0, 2.0, 5.0, 2.0, 7.0, 9.0, 10.0, 6.0, 5.0, 5.0, 4.0, 1.0, 2.0, 3.0, 8.0, 7.0, 7.0, 1.0, 2.0, 6.0, 7.0, 3.0, 8.0, 2.0, 4.0, 5.0, 7.0, 1.0, 2.0, 7.0, 7.0, 5.0, 3.0, 7.0, 6.0, 6.0, 10.0, 3.0, 1.0, 3.0, 4.0, 10.0, 8.0, 10.0, 3.0, 5.0, 4.0, 5.0, 1.0, 3.0, 4.0, 9.0, 2.0, 7.0, 10.0, 6.0, 8.0, 5.0, 8.0, 7.0, 1.0, 3.0, 2.0, 5.0, 9.0, 2.0, 2.0, 1.0, 1.0, 8.0, 4.0, 1.0, 2.0, 2.0, 4.0, 8.0, 9.0, 4.0, 3.0, 3.0, 1.0, 1.0, 3.0, 10.0, 5.0, 5.0, 3.0, 5.0, 9.0, 7.0, 8.0, 8.0, 7.0, 8.0, 1.0, 4.0, 8.0, 5.0, 3.0, 8.0, 4.0, 6.0, 2.0, 9.0, 10.0, 3.0, 9.0, 9.0, 5.0, 2.0, 3.0, 2.0]
global b_x = 5
global d_y = [9.0, 5.0, 9.0, 7.0, 3.0, 4.0, 2.0, 7.0, 9.0, 5.0, 3.0, 10.0, 2.0, 4.0, 4.0, 5.0, 6.0, 4.0, 3.0, 1.0, 2.0, 3.0, 6.0, 10.0, 8.0, 3.0, 5.0, 5.0, 6.0, 4.0, 8.0, 10.0, 6.0, 4.0, 3.0, 4.0, 10.0, 1.0, 7.0, 3.0, 9.0, 3.0, 5.0, 4.0, 5.0, 6.0, 10.0, 1.0, 3.0, 10.0, 10.0, 9.0, 4.0, 3.0, 6.0, 5.0, 8.0, 2.0, 4.0, 4.0, 7.0, 9.0, 3.0, 6.0, 4.0, 9.0, 8.0, 5.0, 4.0, 9.0, 5.0, 4.0, 4.0, 9.0, 7.0, 6.0, 9.0, 2.0, 2.0, 4.0, 2.0, 4.0, 6.0, 3.0, 10.0, 9.0, 1.0, 4.0, 6.0, 10.0, 6.0, 3.0, 3.0, 3.0, 3.0, 6.0, 10.0, 1.0, 3.0, 9.0, 8.0, 10.0, 5.0, 2.0, 8.0, 2.0, 3.0, 5.0, 10.0, 1.0, 8.0, 6.0, 1.0, 7.0, 10.0, 9.0, 7.0, 10.0, 10.0, 6.0, 4.0, 4.0, 2.0, 6.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.437, 0.544, 0.465, 0.525, 0.474, 0.475, 0.427, 0.597, 0.467, 0.589, 0.473, 0.594, 0.596, 0.421, 0.532, 0.553, 0.482, 0.404, 0.424, 0.447, 0.568, 0.443, 0.511, 0.439, 0.55, 0.424, 0.42, 0.556, 0.582, 0.426, 0.505, 0.592, 0.405, 0.544, 0.576, 0.434, 0.409, 0.415, 0.575, 0.48, 0.428, 0.425, 0.504, 0.524, 0.533, 0.407, 0.419, 0.507, 0.451, 0.422, 0.583, 0.547, 0.533, 0.532, 0.495, 0.589, 0.536, 0.518, 0.516, 0.438, 0.464, 0.543, 0.533, 0.531, 0.436, 0.565, 0.524, 0.424, 0.47, 0.419, 0.584, 0.405, 0.42, 0.555, 0.425, 0.598, 0.413, 0.535, 0.519, 0.563, 0.402, 0.473, 0.49, 0.497, 0.411, 0.51, 0.572, 0.512, 0.526, 0.57, 0.596, 0.585, 0.416, 0.535, 0.438, 0.536, 0.529, 0.48, 0.541, 0.578, 0.578, 0.426, 0.433, 0.533, 0.579, 0.561, 0.55, 0.549, 0.506, 0.558, 0.418, 0.438, 0.44, 0.404, 0.516, 0.405, 0.559, 0.554, 0.548, 0.403, 0.507, 0.517, 0.487, 0.433]
global origin = 1
global destination = 35