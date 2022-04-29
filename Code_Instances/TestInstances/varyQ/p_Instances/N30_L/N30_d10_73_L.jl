global arcs = [1 3; 1 5; 1 6; 1 8; 1 10; 2 16; 2 22; 2 23; 2 27; 3 26; 4 18; 4 24; 5 14; 5 15; 5 22; 6 10; 6 18; 6 30; 7 6; 7 11; 7 20; 8 9; 8 12; 8 18; 8 20; 9 4; 9 10; 9 19; 9 21; 9 23; 9 29; 9 30; 10 8; 10 20; 10 29; 11 8; 11 12; 12 8; 12 13; 12 23; 12 29; 13 8; 13 18; 13 19; 13 21; 14 7; 14 13; 14 18; 14 28; 15 7; 16 8; 16 11; 16 20; 17 2; 17 12; 17 21; 17 28; 18 22; 19 5; 19 9; 19 10; 19 14; 19 23; 20 7; 20 25; 21 8; 21 9; 21 18; 21 25; 21 29; 22 25; 23 4; 23 6; 23 8; 23 24; 23 26; 23 28; 24 17; 25 8; 25 10; 25 20; 25 21; 25 26; 25 27; 25 30; 26 5; 26 7; 26 29; 26 30; 27 2; 27 23; 27 26; 28 2; 28 11; 28 16; 28 17; 28 18; 29 3; 29 12; 29 16; 29 19; 29 27]
global d_x = [5.0, 5.0, 10.0, 4.0, 5.0, 2.0, 9.0, 8.0, 3.0, 9.0, 6.0, 7.0, 8.0, 4.0, 9.0, 3.0, 8.0, 4.0, 5.0, 10.0, 9.0, 1.0, 5.0, 6.0, 9.0, 9.0, 4.0, 5.0, 9.0, 4.0, 8.0, 1.0, 9.0, 1.0, 3.0, 5.0, 8.0, 1.0, 2.0, 3.0, 6.0, 1.0, 4.0, 5.0, 2.0, 7.0, 9.0, 8.0, 9.0, 3.0, 2.0, 1.0, 7.0, 3.0, 7.0, 9.0, 3.0, 10.0, 2.0, 10.0, 5.0, 6.0, 4.0, 8.0, 5.0, 10.0, 2.0, 5.0, 10.0, 6.0, 1.0, 9.0, 2.0, 9.0, 9.0, 5.0, 1.0, 8.0, 3.0, 2.0, 7.0, 2.0, 10.0, 10.0, 8.0, 4.0, 1.0, 7.0, 4.0, 7.0, 9.0, 3.0, 10.0, 8.0, 2.0, 8.0, 7.0, 6.0, 9.0, 8.0, 4.0, 10.0]
global b_x = 5
global d_y = [2.0, 6.0, 4.0, 6.0, 1.0, 3.0, 5.0, 8.0, 3.0, 3.0, 7.0, 3.0, 1.0, 7.0, 9.0, 2.0, 1.0, 4.0, 3.0, 3.0, 4.0, 6.0, 6.0, 2.0, 9.0, 10.0, 4.0, 3.0, 9.0, 4.0, 6.0, 3.0, 9.0, 2.0, 1.0, 2.0, 4.0, 4.0, 5.0, 3.0, 5.0, 5.0, 7.0, 5.0, 3.0, 9.0, 1.0, 3.0, 6.0, 6.0, 9.0, 4.0, 3.0, 4.0, 2.0, 3.0, 1.0, 1.0, 6.0, 2.0, 4.0, 3.0, 7.0, 6.0, 3.0, 3.0, 6.0, 6.0, 3.0, 6.0, 5.0, 7.0, 5.0, 1.0, 4.0, 8.0, 6.0, 7.0, 7.0, 1.0, 6.0, 7.0, 8.0, 9.0, 9.0, 10.0, 1.0, 9.0, 1.0, 7.0, 6.0, 10.0, 2.0, 9.0, 2.0, 1.0, 5.0, 8.0, 10.0, 9.0, 2.0, 4.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.584, 0.408, 0.569, 0.423, 0.542, 0.454, 0.589, 0.531, 0.55, 0.552, 0.434, 0.411, 0.45, 0.557, 0.529, 0.545, 0.478, 0.441, 0.455, 0.489, 0.461, 0.45, 0.402, 0.437, 0.583, 0.45, 0.553, 0.537, 0.403, 0.535, 0.442, 0.51, 0.538, 0.555, 0.516, 0.423, 0.516, 0.579, 0.464, 0.514, 0.502, 0.535, 0.417, 0.543, 0.52, 0.447, 0.451, 0.504, 0.564, 0.47, 0.536, 0.442, 0.425, 0.504, 0.477, 0.596, 0.439, 0.405, 0.591, 0.41, 0.407, 0.561, 0.415, 0.501, 0.402, 0.435, 0.492, 0.41, 0.459, 0.437, 0.548, 0.433, 0.502, 0.55, 0.405, 0.492, 0.571, 0.403, 0.478, 0.527, 0.522, 0.514, 0.409, 0.434, 0.435, 0.548, 0.49, 0.48, 0.525, 0.435, 0.434, 0.571, 0.572, 0.577, 0.414, 0.492, 0.502, 0.567, 0.553, 0.511, 0.566, 0.504]
global origin = 1
global destination = 30