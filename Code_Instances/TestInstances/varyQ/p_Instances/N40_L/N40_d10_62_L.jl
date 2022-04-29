global arcs = [1 5; 1 14; 1 24; 2 3; 2 15; 2 30; 2 40; 3 36; 4 8; 5 4; 5 8; 5 11; 5 13; 5 18; 5 19; 5 21; 5 24; 6 7; 6 22; 6 35; 7 8; 7 9; 7 14; 7 17; 7 22; 7 24; 7 30; 8 5; 8 25; 9 30; 10 2; 10 17; 10 28; 10 37; 10 40; 11 6; 11 7; 11 8; 11 39; 12 4; 12 16; 12 20; 13 4; 13 17; 13 20; 14 26; 14 37; 15 5; 15 16; 15 38; 16 5; 16 7; 16 19; 16 25; 16 27; 17 8; 17 15; 17 23; 17 27; 17 30; 17 32; 18 2; 18 7; 18 9; 18 13; 18 20; 18 39; 19 6; 19 12; 19 31; 20 23; 21 7; 21 8; 21 9; 21 26; 21 33; 22 14; 22 28; 23 4; 23 18; 23 21; 23 25; 23 39; 23 40; 24 23; 24 34; 25 6; 25 14; 25 23; 26 13; 26 16; 26 29; 26 30; 26 39; 27 4; 28 31; 28 35; 28 40; 29 12; 30 2; 30 11; 30 17; 31 2; 31 8; 31 18; 32 8; 32 31; 32 33; 32 40; 33 19; 33 20; 33 27; 34 15; 34 22; 34 33; 35 11; 35 12; 36 4; 36 8; 36 14; 36 28; 36 30; 36 33; 36 37; 37 23; 37 33; 38 5; 38 13; 38 14; 38 22; 38 25; 38 36; 39 10; 39 12; 39 17; 39 40]
global d_x = [6.0, 4.0, 9.0, 9.0, 5.0, 2.0, 8.0, 3.0, 2.0, 3.0, 5.0, 1.0, 7.0, 10.0, 3.0, 7.0, 10.0, 6.0, 6.0, 4.0, 9.0, 6.0, 10.0, 7.0, 8.0, 9.0, 2.0, 4.0, 2.0, 4.0, 2.0, 8.0, 9.0, 7.0, 9.0, 1.0, 8.0, 1.0, 9.0, 3.0, 6.0, 8.0, 3.0, 6.0, 3.0, 3.0, 3.0, 3.0, 1.0, 4.0, 2.0, 3.0, 2.0, 4.0, 9.0, 1.0, 2.0, 6.0, 7.0, 4.0, 4.0, 4.0, 4.0, 3.0, 8.0, 9.0, 5.0, 1.0, 9.0, 2.0, 2.0, 5.0, 5.0, 3.0, 1.0, 9.0, 8.0, 2.0, 1.0, 2.0, 10.0, 3.0, 5.0, 2.0, 1.0, 1.0, 10.0, 3.0, 10.0, 5.0, 10.0, 6.0, 4.0, 7.0, 9.0, 1.0, 5.0, 7.0, 2.0, 4.0, 7.0, 6.0, 1.0, 9.0, 5.0, 2.0, 4.0, 3.0, 1.0, 2.0, 2.0, 6.0, 4.0, 8.0, 3.0, 8.0, 6.0, 6.0, 6.0, 7.0, 4.0, 8.0, 8.0, 6.0, 3.0, 7.0, 9.0, 6.0, 2.0, 2.0, 2.0, 4.0, 1.0, 4.0, 9.0, 1.0]
global b_x = 5
global d_y = [8.0, 4.0, 1.0, 10.0, 2.0, 8.0, 2.0, 1.0, 1.0, 6.0, 10.0, 6.0, 6.0, 10.0, 5.0, 6.0, 4.0, 3.0, 1.0, 8.0, 10.0, 9.0, 10.0, 5.0, 10.0, 5.0, 1.0, 4.0, 7.0, 6.0, 6.0, 3.0, 7.0, 8.0, 8.0, 3.0, 10.0, 10.0, 10.0, 6.0, 8.0, 4.0, 2.0, 4.0, 3.0, 4.0, 1.0, 1.0, 3.0, 3.0, 6.0, 6.0, 4.0, 7.0, 7.0, 2.0, 9.0, 9.0, 4.0, 7.0, 3.0, 2.0, 9.0, 8.0, 8.0, 6.0, 2.0, 3.0, 5.0, 4.0, 10.0, 8.0, 9.0, 3.0, 10.0, 1.0, 3.0, 9.0, 5.0, 3.0, 9.0, 7.0, 6.0, 2.0, 8.0, 10.0, 10.0, 1.0, 2.0, 7.0, 4.0, 1.0, 8.0, 5.0, 1.0, 2.0, 8.0, 8.0, 5.0, 7.0, 7.0, 1.0, 4.0, 9.0, 3.0, 8.0, 5.0, 2.0, 3.0, 8.0, 1.0, 9.0, 3.0, 9.0, 2.0, 1.0, 7.0, 9.0, 2.0, 8.0, 1.0, 8.0, 2.0, 3.0, 6.0, 10.0, 1.0, 9.0, 8.0, 9.0, 10.0, 4.0, 4.0, 2.0, 4.0, 9.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.489, 0.466, 0.424, 0.505, 0.403, 0.422, 0.406, 0.527, 0.408, 0.583, 0.582, 0.539, 0.403, 0.559, 0.582, 0.428, 0.465, 0.554, 0.534, 0.522, 0.58, 0.519, 0.415, 0.556, 0.573, 0.582, 0.405, 0.542, 0.577, 0.595, 0.473, 0.407, 0.497, 0.5, 0.405, 0.514, 0.432, 0.573, 0.456, 0.502, 0.433, 0.53, 0.552, 0.461, 0.453, 0.53, 0.457, 0.513, 0.51, 0.581, 0.557, 0.56, 0.419, 0.598, 0.454, 0.485, 0.416, 0.435, 0.45, 0.531, 0.495, 0.485, 0.448, 0.42, 0.554, 0.535, 0.445, 0.567, 0.409, 0.461, 0.524, 0.57, 0.445, 0.449, 0.421, 0.45, 0.437, 0.45, 0.561, 0.522, 0.431, 0.54, 0.453, 0.58, 0.505, 0.471, 0.463, 0.522, 0.546, 0.477, 0.411, 0.426, 0.418, 0.525, 0.514, 0.446, 0.493, 0.509, 0.425, 0.426, 0.533, 0.494, 0.407, 0.416, 0.568, 0.577, 0.415, 0.562, 0.593, 0.431, 0.569, 0.425, 0.54, 0.469, 0.553, 0.536, 0.452, 0.565, 0.444, 0.495, 0.499, 0.494, 0.489, 0.497, 0.421, 0.586, 0.449, 0.537, 0.516, 0.478, 0.478, 0.528, 0.44, 0.57, 0.583, 0.561]
global origin = 1
global destination = 40