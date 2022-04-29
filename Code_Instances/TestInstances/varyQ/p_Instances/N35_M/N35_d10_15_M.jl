global arcs = [1 2; 1 4; 1 10; 1 20; 1 35; 2 14; 2 18; 2 22; 2 29; 3 4; 3 11; 3 13; 3 29; 3 32; 4 7; 4 8; 4 15; 4 24; 4 25; 4 29; 4 32; 5 32; 6 10; 6 11; 6 28; 7 3; 7 17; 8 2; 8 9; 8 21; 8 25; 8 34; 9 5; 9 10; 9 25; 10 18; 10 31; 11 8; 11 10; 11 12; 11 25; 11 27; 11 28; 11 32; 12 3; 12 8; 12 11; 12 14; 12 19; 12 21; 13 25; 14 31; 14 35; 15 19; 15 22; 16 2; 16 23; 17 5; 17 13; 17 27; 17 35; 18 6; 18 14; 18 16; 19 23; 20 4; 20 22; 20 25; 20 31; 21 4; 21 10; 21 15; 22 2; 22 18; 22 21; 22 32; 22 33; 22 34; 23 11; 23 12; 23 18; 23 22; 23 27; 24 9; 24 26; 24 28; 25 3; 25 5; 25 14; 25 31; 26 5; 26 10; 26 14; 26 33; 27 18; 27 24; 28 8; 28 11; 28 24; 28 33; 28 34; 29 31; 30 2; 30 12; 30 25; 30 26; 31 4; 31 19; 31 34; 32 2; 32 5; 32 7; 32 14; 32 21; 32 30; 33 32; 34 15; 34 19; 34 25]
global d_x = [6.0, 8.0, 2.0, 4.0, 10.0, 4.0, 2.0, 10.0, 3.0, 10.0, 3.0, 8.0, 7.0, 2.0, 8.0, 1.0, 6.0, 2.0, 3.0, 8.0, 7.0, 6.0, 3.0, 4.0, 8.0, 8.0, 5.0, 7.0, 2.0, 5.0, 10.0, 1.0, 9.0, 6.0, 4.0, 6.0, 3.0, 6.0, 7.0, 6.0, 3.0, 9.0, 7.0, 7.0, 7.0, 9.0, 6.0, 5.0, 10.0, 2.0, 5.0, 10.0, 10.0, 9.0, 8.0, 8.0, 6.0, 5.0, 6.0, 9.0, 2.0, 5.0, 8.0, 5.0, 2.0, 1.0, 8.0, 8.0, 2.0, 5.0, 4.0, 4.0, 10.0, 2.0, 7.0, 1.0, 9.0, 9.0, 5.0, 6.0, 10.0, 9.0, 10.0, 3.0, 5.0, 1.0, 8.0, 7.0, 3.0, 1.0, 5.0, 10.0, 6.0, 2.0, 8.0, 3.0, 1.0, 1.0, 8.0, 10.0, 9.0, 2.0, 10.0, 1.0, 7.0, 2.0, 1.0, 10.0, 10.0, 2.0, 2.0, 9.0, 1.0, 10.0, 9.0, 3.0, 9.0, 10.0, 4.0]
global b_x = 5
global d_y = [9.0, 1.0, 2.0, 4.0, 10.0, 5.0, 8.0, 1.0, 1.0, 7.0, 1.0, 8.0, 10.0, 10.0, 3.0, 5.0, 3.0, 7.0, 1.0, 8.0, 10.0, 4.0, 3.0, 10.0, 3.0, 9.0, 7.0, 8.0, 6.0, 10.0, 10.0, 2.0, 5.0, 6.0, 7.0, 5.0, 8.0, 1.0, 6.0, 8.0, 6.0, 6.0, 8.0, 2.0, 5.0, 8.0, 9.0, 5.0, 8.0, 4.0, 10.0, 1.0, 7.0, 5.0, 2.0, 8.0, 8.0, 9.0, 4.0, 6.0, 2.0, 9.0, 10.0, 5.0, 2.0, 1.0, 10.0, 6.0, 8.0, 3.0, 9.0, 9.0, 3.0, 6.0, 6.0, 3.0, 6.0, 5.0, 3.0, 6.0, 5.0, 3.0, 1.0, 10.0, 4.0, 3.0, 10.0, 7.0, 5.0, 2.0, 10.0, 10.0, 6.0, 9.0, 3.0, 6.0, 9.0, 9.0, 9.0, 1.0, 3.0, 7.0, 6.0, 5.0, 2.0, 2.0, 7.0, 3.0, 2.0, 1.0, 10.0, 1.0, 4.0, 7.0, 10.0, 5.0, 5.0, 4.0, 5.0]
global b_y = 10
global p = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
global q = [0.684, 0.736, 0.674, 0.748, 0.659, 0.717, 0.661, 0.76, 0.669, 0.761, 0.661, 0.619, 0.623, 0.616, 0.795, 0.673, 0.639, 0.736, 0.655, 0.684, 0.653, 0.766, 0.675, 0.762, 0.745, 0.769, 0.792, 0.689, 0.692, 0.63, 0.762, 0.606, 0.711, 0.769, 0.605, 0.753, 0.611, 0.697, 0.73, 0.657, 0.628, 0.678, 0.724, 0.684, 0.729, 0.601, 0.638, 0.689, 0.707, 0.648, 0.784, 0.616, 0.717, 0.618, 0.734, 0.688, 0.789, 0.755, 0.644, 0.759, 0.797, 0.618, 0.781, 0.649, 0.66, 0.703, 0.774, 0.77, 0.71, 0.728, 0.679, 0.736, 0.787, 0.608, 0.681, 0.613, 0.608, 0.788, 0.788, 0.76, 0.745, 0.621, 0.743, 0.63, 0.674, 0.691, 0.609, 0.78, 0.791, 0.76, 0.774, 0.749, 0.68, 0.755, 0.609, 0.624, 0.644, 0.686, 0.615, 0.613, 0.691, 0.64, 0.78, 0.704, 0.754, 0.709, 0.753, 0.768, 0.736, 0.667, 0.708, 0.687, 0.686, 0.737, 0.612, 0.705, 0.626, 0.763, 0.722]
global origin = 1
global destination = 35