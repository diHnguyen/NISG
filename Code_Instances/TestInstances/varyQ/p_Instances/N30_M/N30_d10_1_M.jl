global arcs = [1 12; 1 14; 2 11; 2 18; 2 20; 3 6; 3 10; 3 20; 3 23; 3 25; 3 29; 3 30; 4 16; 4 19; 5 21; 5 28; 6 9; 6 11; 6 14; 6 22; 7 2; 7 8; 7 29; 8 17; 8 23; 8 30; 9 13; 10 9; 10 15; 11 7; 11 22; 11 23; 11 24; 11 30; 12 2; 12 14; 12 24; 12 27; 13 3; 13 6; 14 10; 14 11; 14 15; 14 18; 14 19; 14 22; 14 25; 14 27; 15 2; 15 6; 15 9; 15 24; 16 13; 17 22; 17 23; 17 28; 18 16; 19 5; 19 10; 20 27; 21 6; 21 15; 22 2; 22 4; 22 27; 23 12; 24 4; 24 16; 24 26; 25 5; 25 21; 26 3; 26 4; 26 16; 26 28; 26 30; 27 10; 27 16; 27 24; 27 25; 27 29; 28 4; 28 5; 28 7; 28 8; 28 25; 29 8; 29 18; 29 27; 29 28]
global d_x = [5.0, 5.0, 10.0, 2.0, 4.0, 3.0, 3.0, 4.0, 2.0, 1.0, 9.0, 6.0, 1.0, 4.0, 4.0, 9.0, 2.0, 1.0, 5.0, 1.0, 7.0, 5.0, 5.0, 6.0, 3.0, 3.0, 4.0, 10.0, 3.0, 3.0, 7.0, 4.0, 4.0, 3.0, 10.0, 6.0, 10.0, 5.0, 9.0, 2.0, 6.0, 7.0, 6.0, 6.0, 4.0, 5.0, 8.0, 6.0, 1.0, 1.0, 9.0, 10.0, 5.0, 3.0, 3.0, 5.0, 6.0, 3.0, 9.0, 10.0, 7.0, 1.0, 4.0, 3.0, 4.0, 9.0, 7.0, 9.0, 9.0, 2.0, 7.0, 6.0, 6.0, 1.0, 5.0, 9.0, 1.0, 10.0, 7.0, 1.0, 5.0, 7.0, 10.0, 9.0, 6.0, 7.0, 1.0, 9.0, 2.0, 5.0]
global b_x = 5
global d_y = [4.0, 10.0, 8.0, 6.0, 8.0, 2.0, 8.0, 7.0, 8.0, 9.0, 9.0, 5.0, 3.0, 4.0, 10.0, 9.0, 7.0, 7.0, 3.0, 10.0, 3.0, 3.0, 1.0, 8.0, 3.0, 4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 4.0, 3.0, 5.0, 10.0, 1.0, 8.0, 6.0, 3.0, 1.0, 1.0, 4.0, 6.0, 8.0, 7.0, 4.0, 3.0, 1.0, 10.0, 7.0, 3.0, 5.0, 6.0, 9.0, 10.0, 7.0, 7.0, 9.0, 8.0, 2.0, 8.0, 3.0, 3.0, 6.0, 9.0, 1.0, 7.0, 1.0, 2.0, 10.0, 10.0, 4.0, 8.0, 3.0, 7.0, 6.0, 7.0, 5.0, 7.0, 9.0, 3.0, 10.0, 10.0, 5.0, 5.0, 10.0, 7.0, 6.0, 3.0, 3.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.693, 0.714, 0.663, 0.741, 0.798, 0.749, 0.67, 0.687, 0.666, 0.629, 0.609, 0.624, 0.636, 0.623, 0.656, 0.657, 0.775, 0.719, 0.792, 0.679, 0.631, 0.783, 0.679, 0.622, 0.763, 0.672, 0.726, 0.6, 0.658, 0.656, 0.767, 0.647, 0.766, 0.621, 0.741, 0.781, 0.751, 0.734, 0.71, 0.69, 0.783, 0.771, 0.727, 0.733, 0.632, 0.74, 0.651, 0.691, 0.717, 0.796, 0.743, 0.721, 0.776, 0.624, 0.704, 0.729, 0.619, 0.652, 0.639, 0.654, 0.793, 0.775, 0.724, 0.61, 0.688, 0.681, 0.764, 0.736, 0.731, 0.693, 0.635, 0.644, 0.611, 0.796, 0.792, 0.781, 0.629, 0.737, 0.77, 0.62, 0.605, 0.671, 0.728, 0.659, 0.69, 0.665, 0.752, 0.668, 0.691, 0.669]
global origin = 1
global destination = 30