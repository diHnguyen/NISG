global arcs = [1 10; 1 17; 1 30; 2 3; 2 23; 2 25; 3 11; 3 12; 3 16; 3 26; 4 2; 4 3; 4 12; 4 17; 4 24; 5 4; 5 22; 5 26; 6 17; 6 25; 6 29; 7 9; 7 26; 7 28; 8 9; 8 11; 8 16; 8 25; 9 6; 9 16; 9 26; 10 13; 10 22; 11 9; 11 19; 11 25; 12 3; 12 18; 12 23; 13 3; 13 5; 13 9; 14 12; 14 19; 15 25; 15 29; 16 5; 16 24; 17 6; 17 13; 17 20; 17 21; 17 24; 17 25; 17 26; 17 27; 18 5; 18 8; 18 12; 18 13; 18 16; 18 23; 18 25; 19 4; 19 17; 19 30; 20 9; 20 22; 20 26; 21 9; 21 12; 21 25; 22 10; 23 5; 23 15; 24 7; 24 20; 25 5; 26 4; 26 27; 26 30; 27 12; 28 2; 28 14; 28 15; 28 21; 28 24; 28 26; 28 27; 29 17]
global d_x = [3.0, 3.0, 8.0, 2.0, 10.0, 3.0, 1.0, 8.0, 9.0, 2.0, 10.0, 6.0, 8.0, 10.0, 5.0, 3.0, 5.0, 8.0, 7.0, 4.0, 3.0, 2.0, 6.0, 10.0, 9.0, 2.0, 9.0, 2.0, 3.0, 1.0, 1.0, 2.0, 5.0, 6.0, 1.0, 3.0, 5.0, 5.0, 2.0, 9.0, 1.0, 7.0, 10.0, 6.0, 2.0, 9.0, 10.0, 2.0, 10.0, 2.0, 9.0, 2.0, 8.0, 2.0, 1.0, 10.0, 4.0, 3.0, 6.0, 4.0, 7.0, 3.0, 8.0, 2.0, 6.0, 3.0, 2.0, 5.0, 7.0, 2.0, 6.0, 5.0, 8.0, 5.0, 5.0, 5.0, 5.0, 9.0, 7.0, 5.0, 7.0, 1.0, 2.0, 8.0, 10.0, 5.0, 3.0, 10.0, 4.0, 7.0]
global b_x = 5
global d_y = [2.0, 10.0, 4.0, 3.0, 5.0, 3.0, 5.0, 3.0, 2.0, 1.0, 2.0, 7.0, 8.0, 5.0, 5.0, 3.0, 6.0, 4.0, 4.0, 5.0, 3.0, 7.0, 9.0, 5.0, 6.0, 7.0, 10.0, 6.0, 5.0, 3.0, 3.0, 6.0, 5.0, 8.0, 9.0, 6.0, 9.0, 8.0, 6.0, 9.0, 4.0, 1.0, 7.0, 1.0, 10.0, 5.0, 10.0, 8.0, 6.0, 9.0, 4.0, 1.0, 3.0, 9.0, 1.0, 8.0, 9.0, 9.0, 8.0, 1.0, 1.0, 2.0, 8.0, 2.0, 4.0, 6.0, 8.0, 5.0, 8.0, 1.0, 7.0, 5.0, 1.0, 6.0, 2.0, 8.0, 3.0, 9.0, 10.0, 9.0, 2.0, 5.0, 4.0, 5.0, 8.0, 10.0, 3.0, 2.0, 5.0, 9.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.767, 0.62, 0.67, 0.628, 0.722, 0.612, 0.714, 0.61, 0.666, 0.618, 0.655, 0.72, 0.615, 0.748, 0.638, 0.787, 0.61, 0.715, 0.655, 0.777, 0.622, 0.717, 0.629, 0.628, 0.766, 0.631, 0.723, 0.758, 0.613, 0.78, 0.761, 0.688, 0.745, 0.74, 0.647, 0.671, 0.679, 0.601, 0.796, 0.636, 0.771, 0.605, 0.606, 0.731, 0.752, 0.675, 0.619, 0.626, 0.701, 0.682, 0.779, 0.687, 0.702, 0.766, 0.796, 0.683, 0.614, 0.759, 0.612, 0.701, 0.673, 0.623, 0.648, 0.729, 0.625, 0.701, 0.603, 0.797, 0.671, 0.734, 0.664, 0.761, 0.782, 0.659, 0.684, 0.69, 0.678, 0.662, 0.793, 0.611, 0.65, 0.702, 0.699, 0.66, 0.621, 0.792, 0.737, 0.626, 0.665, 0.726]
global origin = 1
global destination = 30